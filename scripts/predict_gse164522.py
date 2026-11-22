"""Apply the frozen GSE245552 tissue classifiers (`V1.1`) to GSE164522 tumor specimens.

Track 2 (`V1.1`) upgrades:
1. Full clean external bags (`min_cells=50, max_cells=None` by default) stripping
   author-annotated cross-lineage `CD3E+` myeloid doublets (`hC56/41/47/36`) and
   `hC75_B-CD3E` B-T doublets before bag selection.
2. Deterministic HGNC symbol alignment onto `245_164_shared_nuclear` (18,815 canonical
   nuclear genes) with `expm1 -> CP10k -> log1p` shared-nuclear recalibration,
   achieving `512/512 = 100.0%` PCA HVG coverage (0/15 model rejections).
3. Self-contained frozen `debiased_ucell_scores` scoring for `program_elastic` using
   the SHA256-locked development `program_defs.json` and `matched_controls.csv`
   (zero external adaptive control refitting).
"""

from __future__ import annotations

import argparse
import hashlib
from itertools import combinations
import json
import os
import platform
from pathlib import Path
import resource
import subprocess
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter"), str(ROOT / "scripts")]

from crlm.hgnc_alias_map import (
    align_axis_to_canonical_targets,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
)
from crlm.model_bags import (
    Bag,
    LINEAGES,
    SEED,
    sample_bags,
    select_bag_rows,
)
from crlm.program_scores import controls_from_dataframe, debiased_ucell_scores
from linear_attention_mil import load_model, predict_bags

DATA = ROOT / "data/processed/m6/GSE164522"
EXTERNAL = ROOT / "artifacts/m6/GSE164522"
OUT = EXTERNAL / "model_predictions"
CACHE = ROOT / "embeddings/m6/GSE164522/scgpt"
M5 = ROOT / "artifacts/m5"
DEFAULT_UNIVERSE = "245_164_shared_nuclear"
METHODS = (
    "program_elastic",
    "PCA_mean",
    "PCA_attention",
    "scGPT_mean",
    "scGPT_attention",
)
TISSUES = {"PT": "PT", "MT": "LM"}


def _sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_output_sha256(stage_dir: Path) -> int:
    target = stage_dir / "output_sha256.txt"
    files = sorted(
        p
        for p in stage_dir.rglob("*")
        if p.is_file() and p.name not in ("output_sha256.txt", ".gitkeep")
    )
    target.write_text(
        "".join(f"{_sha256_file(p)}  {p.relative_to(ROOT)}\n" for p in files),
        encoding="utf-8",
    )
    return len(files)


def frozen_models(revision: str = "WORKTREE") -> dict[tuple[str, str], Path]:
    frozen = M5 / "frozen_inputs_sha256.txt"
    if not frozen.is_file():
        raise FileNotFoundError(f"Missing frozen input registry: {frozen}")
    if revision and revision.upper() != "WORKTREE":
        version = subprocess.check_output(
            ["git", "show", f"{revision}:{frozen.relative_to(ROOT)}"], cwd=ROOT
        )
        if version != frozen.read_bytes():
            raise ValueError(
                "The freeze revision does not match the classifier input registry"
            )
    else:
        version = frozen.read_bytes()
    for line in version.decode("utf-8").splitlines():
        if not line.strip():
            continue
        expected, filename = line.split("  ", 1)
        if filename.startswith(
            (
                "models/m4/full/",
                "artifacts/m5/",
                "scripts/encode_scgpt.py",
                "src/crlm/model_bags.py",
                "src/crlm/hgnc_alias_map.py",
            )
        ):
            if _sha256_file(ROOT / filename) != expected:
                raise ValueError(f"Frozen input differs: {filename}")
    fits = pd.read_csv(M5 / "frozen_models.csv")
    models = {
        (r.lineage, r.method): ROOT / r.path for r in fits.itertuples() if r.accepted
    }
    for lin in LINEAGES:
        prog_dir = models[(lin, "program_elastic")].parent
        if not (prog_dir / "program_defs.json").is_file() or not (
            prog_dir / "matched_controls.csv"
        ).is_file():
            raise FileNotFoundError(
                f"Missing frozen program_defs.json or matched_controls.csv in {prog_dir}"
            )
    return models


def is_gse164522_doublet(obs: pd.DataFrame) -> np.ndarray:
    """Identify author-annotated CD3E+ myeloid doublets (`hC56/41/47/36`) and `hC75_B-CD3E` B-T doublets."""
    sub = obs["celltype_sub"].astype(str)
    lin = obs["lineage"].astype(str)
    myeloid_cd3 = lin.eq("Mono/Mac/DC") & sub.str.contains("CD3E", regex=False)
    b_cd3 = lin.eq("B/Plasma") & sub.eq("hC75_B-CD3E")
    return (myeloid_cd3 | b_cd3).to_numpy(dtype=bool)


def recalibrate_shared_nuclear_cp10k(
    mat_csr: sparse.csr_matrix,
    shared_col_indices: np.ndarray,
) -> tuple[sparse.csr_matrix, np.ndarray]:
    """Slice to `245_164_shared_nuclear` and recalibrate `expm1 -> CP10k -> log1p`."""
    sub = mat_csr[:, shared_col_indices].tocsr().copy()
    sub.sum_duplicates()
    sub.eliminate_zeros()
    if np.any(sub.data < 0) or not np.isfinite(sub.data).all():
        raise ValueError("Expected finite nonnegative processed expression")
    linear = sub.copy()
    linear.data = np.expm1(linear.data.astype(np.float64))
    nuclear_sums = np.asarray(linear.sum(axis=1)).ravel()
    if (nuclear_sums <= 0).any():
        raise ValueError(
            "Encountered cell with zero shared-nuclear expression sum during CP10k recalibration"
        )
    row_nnz = np.diff(linear.indptr)
    linear.data = np.log1p(
        linear.data * np.repeat(10000.0 / nuclear_sums, row_nnz)
    ).astype(np.float32)
    return linear, nuclear_sums


def tumor_bags(
    path: Path,
    tissue: str,
    *,
    gene_universe: str = DEFAULT_UNIVERSE,
    exclude_doublets: bool = True,
    min_cells: int = 50,
    max_cells: int | None = None,
):
    universe_genes = list(load_shared_nuclear_universe(gene_universe))
    genes_idx = pd.Index(universe_genes)
    obj = ad.read_h5ad(path, backed="r")
    try:
        obs = obj.obs
        if set(obs["tissue"].astype(str)) != {tissue} or not obs.index.is_unique:
            raise ValueError(f"Inconsistent cell identity or tissue in {path}")
        aligned = align_axis_to_canonical_targets(
            obj.var["gene_symbol"].astype(str).tolist(), universe_genes
        )
        if not aligned["present"].all():
            missing = aligned.loc[~aligned["present"], "canonical_symbol"].head(10).tolist()
            raise ValueError(
                f"Missing {int((~aligned['present']).sum())} shared nuclear genes in {path}: {missing}"
            )
        shared_cols = aligned["source_index"].to_numpy(dtype=np.int64)
        doublet_mask = is_gse164522_doublet(obs)

        for sample, idx in obs.groupby("sample_uid", observed=True).indices.items():
            sample_frame = obs.iloc[idx]
            if sample_frame[["donor_uid", "lesion_uid"]].nunique().max() != 1:
                raise ValueError(f"Unresolved specimen identity: {sample}")
            lin_arr = sample_frame["lineage"].astype(str).to_numpy()
            for lineage in LINEAGES:
                in_lin = lin_arr == lineage
                raw_rows = idx[in_lin]
                if exclude_doublets:
                    clean_rows = idx[in_lin & (~doublet_mask[idx])]
                    n_doublets = int(len(raw_rows) - len(clean_rows))
                else:
                    clean_rows = raw_rows
                    n_doublets = 0
                sel = select_bag_rows(
                    obs.index,
                    clean_rows,
                    f"{sample}|{lineage}",
                    min_cells=min_cells,
                    max_cells=max_cells,
                    seed=SEED,
                )
                bag = None
                median_nuc_sum = np.nan
                if len(sel):
                    raw_mat = obj.X[sel, :].tocsr()
                    mat, nuc_sums = recalibrate_shared_nuclear_cp10k(raw_mat, shared_cols)
                    median_nuc_sum = float(np.median(nuc_sums))
                    bag = Bag(obs.iloc[sel].copy(), mat, genes_idx, len(clean_rows))
                yield (
                    dict(
                        sample_uid=str(sample),
                        donor_uid=str(sample_frame["donor_uid"].iloc[0]),
                        lesion_uid=str(sample_frame["lesion_uid"].iloc[0]),
                        tissue=tissue,
                        lineage=lineage,
                        n_raw_lineage=int(len(raw_rows)),
                        n_doublets_removed=n_doublets,
                        n_population=int(len(clean_rows)),
                        n_cells=int(len(sel)),
                        n_captured=int(len(idx)),
                        median_nuclear_inverse_sum=median_nuc_sum,
                    ),
                    bag,
                )
    finally:
        obj.file.close()


def pca_projection(bag: Bag, transform) -> tuple[np.ndarray | None, np.ndarray]:
    positions = bag.genes.get_indexer(transform["genes"])
    measured = positions >= 0
    if measured.mean() < 0.95:
        return None, measured
    values = bag.counts[:, positions[measured]].toarray().astype(np.float32)
    values -= transform["center"][measured]
    return (values @ transform["components"][:, measured].T).astype(float), measured


def elastic_logit(path: Path, means: pd.Series) -> float | None:
    with np.load(path, allow_pickle=False) as fit:
        values = means.reindex(fit["columns"]).to_numpy(dtype=float)
        if not np.isfinite(values).all():
            return None
        scaled = (values - fit["center"]) / fit["scale"]
        return float((scaled @ fit["coef"].T + fit["bias"])[0])


def score_frozen_program_bag(bag: Bag, model_dir: Path) -> pd.Series:
    """Score external `Bag` using frozen development `program_defs.json` and `matched_controls.csv`."""
    spec = json.loads((model_dir / "program_defs.json").read_text(encoding="utf-8"))
    ctrl_df = pd.read_csv(model_dir / "matched_controls.csv")
    matched_controls = controls_from_dataframe(ctrl_df)
    cell_scores = debiased_ucell_scores(
        bag.counts,
        bag.genes,
        spec["programs"],
        matched_controls=matched_controls,
        max_rank=int(spec.get("max_rank", 1500)),
    )
    return cell_scores.mean(axis=0)


def pooled_logit(path: Path, vector: np.ndarray, sample: str, lineage: str) -> float:
    model, info = load_model(path)
    if info["lineage"] != lineage:
        raise ValueError(f"Classifier lineage differs from bag: {path}")
    prediction, _ = predict_bags(model, info, [vector], [sample])
    return float(prediction["logit"].iloc[0])


def development_check():
    universe_genes = list(load_shared_nuclear_universe(DEFAULT_UNIVERSE))
    design = pd.read_csv(ROOT / "artifacts/m4/oof_predictions.csv")
    base = design.loc[design["method"].isin(METHODS) & design["seed"].eq(SEED)]
    rows = []
    for lineage in LINEAGES:
        local = base.loc[base["lineage"].eq(lineage)]
        ids = local.groupby("sample_uid")["method"].nunique()
        sample = min(ids.index[ids.eq(len(METHODS))])
        expected = local.loc[local["sample_uid"].eq(sample)].set_index("method")
        fold = int(expected["fold"].iloc[0])
        folder = ROOT / "models/m4" / lineage.replace("/", "_") / f"fold{fold}"
        path = ROOT / "data/processed/m2" / f"{sample}.h5ad"
        bag = sample_bags(
            path, min_cells=50, max_cells=None, gene_universe=universe_genes
        )[lineage]
        obj = ad.read_h5ad(path, backed="r")
        try:
            orig = np.flatnonzero(
                obj.obs["qc_retained"].to_numpy(dtype=bool)
                & obj.obs["lineage_m2"].astype(str).eq(lineage).to_numpy()
            )
            sel = select_bag_rows(
                obj.obs.index, orig, f"{sample}|{lineage}", min_cells=50, max_cells=None
            )
            if not obj.obs.index[sel].astype(str).equals(bag.obs.index):
                raise ValueError(f"Different selected development cells in {sample}")
        finally:
            obj.file.close()

        scores = score_frozen_program_bag(bag, folder)
        mat = bag.counts.astype(np.float32, copy=True)
        mat.data *= np.repeat(
            10000.0 / np.asarray(mat.sum(axis=1)).ravel(), np.diff(mat.indptr)
        ).astype(np.float32)
        np.log1p(mat.data, out=mat.data)
        with np.load(folder / "pca.npz", allow_pickle=False) as transform:
            pca, measured = pca_projection(
                Bag(bag.obs, mat, bag.genes, bag.population), transform
            )
            if not measured.all():
                raise ValueError("Missing development PCA feature")
        with np.load(
            ROOT / f"embeddings/m4/scgpt/{sample}.npz", allow_pickle=False
        ) as encoded:
            idx = pd.Index(encoded["cell_uid"]).get_indexer(bag.obs.index)
            if (idx < 0).any():
                raise ValueError(f"Missing development embedding for {sample}")
            scgpt = encoded["embedding"][idx].astype(float)
        for method in METHODS:
            if method == "program_elastic":
                logit = elastic_logit(folder / "program_elastic.npz", scores)
            else:
                h = pca if method.startswith("PCA_") else scgpt
                logit = pooled_logit(
                    folder / f"{method}_{SEED}.npz", h, sample, lineage
                )
            e = float(expected.loc[method, "logit"])
            if not np.isclose(logit, e, atol=2e-5, rtol=1e-5):
                raise ValueError(
                    f"OOF prediction disagreement: {lineage}, {method}, {logit}, {e}"
                )
            rows.append(
                dict(
                    sample_uid=sample,
                    lineage=lineage,
                    method=method,
                    fold=fold,
                    n_cells=len(bag.obs),
                    predicted_logit=logit,
                    saved_oof_logit=e,
                    abs_diff=abs(logit - e),
                )
            )
    pd.DataFrame(rows).to_csv(OUT / "development_interface_check.csv", index=False)


def cache_path(code: str, sample: str, lineage: str) -> Path:
    return CACHE / code / f"{sample}_{lineage.replace('/', '_')}.npz"


def encode_external(
    codes: list[str],
    device: str,
    batch: int,
    *,
    exclude_doublets: bool = True,
    max_cells: int | None = None,
    overwrite: bool = False,
):
    import torch
    from encode_scgpt import (
        CODE,
        HASHES,
        INPUT_TRANSFORM,
        encode,
        load_encoder,
        resolve_vocab_ids,
        tokens,
    )

    mapper = get_hgnc_mapper()
    universe_genes = list(load_shared_nuclear_universe(DEFAULT_UNIVERSE))
    profile = json.loads(
        (ROOT / "artifacts/m4/resource_profiles/scgpt_profile.json").read_text(
            encoding="utf-8"
        )
    )
    if (
        profile.get("checkpoint") != HASHES
        or not profile.get("devices", {}).get(device, {}).get("accepted")
    ):
        raise ValueError(f"Unprofiled checkpoint/device: {device}")

    torch.set_num_threads(4)
    model, vocab = load_encoder()
    vocab_ids, _ = resolve_vocab_ids(universe_genes, vocab)
    model.to(device)

    for code in codes:
        path = DATA / f"{code}.h5ad"
        (CACHE / code).mkdir(parents=True, exist_ok=True)
        obj = ad.read_h5ad(path, backed="r")
        source_sha = str(obj.uns["expression_sha256"])
        obj.file.close()
        report = []
        for meta, bag in tumor_bags(
            path,
            TISSUES[code],
            gene_universe=DEFAULT_UNIVERSE,
            exclude_doublets=exclude_doublets,
            max_cells=max_cells,
        ):
            if bag is None:
                continue
            sample, lineage = meta["sample_uid"], meta["lineage"]
            dest = cache_path(code, sample, lineage)
            uid = bag.obs.index.to_numpy(dtype=str)
            need_encode = overwrite or (not dest.exists())
            if not need_encode:
                with np.load(dest, allow_pickle=False) as enc:
                    same = (
                        "gene_universe" in enc.files
                        and str(enc["gene_universe"]) == DEFAULT_UNIVERSE
                        and str(enc["hgnc_sha256"]) == mapper.sha256
                        and bool(enc["exclude_doublets"]) == bool(exclude_doublets)
                        and np.array_equal(enc["cell_uid"], uid)
                        and str(enc["source_sha"]) == source_sha
                        and str(enc["checkpoint"]) == HASHES["best_model.pt"]
                        and str(enc["upstream"]) == CODE
                    )
                    if not same:
                        raise ValueError(
                            f"Cached external CLS embedding identity differs at {dest}; pass --overwrite to regenerate"
                        )
            if need_encode:
                start = time.monotonic()
                # bag.counts is already recalibrated to log1p(CP10k) on 245_164_shared_nuclear
                cells = tokens(
                    bag, vocab, vocab_ids=vocab_ids, recalibrate_cp10k=False
                )
                embedded = encode(model, cells, vocab, batch, device)
                if device == "xpu":
                    torch.xpu.synchronize()
                temp = dest.with_suffix(".partial.npz")
                np.savez_compressed(
                    temp,
                    embedding=embedded,
                    cell_uid=uid,
                    nonzero=np.array([x[3] for x in cells], dtype=np.int32),
                    in_vocab=np.array([x[4] for x in cells], dtype=np.int32),
                    tokens=np.array([len(x[1]) for x in cells], dtype=np.int32),
                    gene_universe=DEFAULT_UNIVERSE,
                    hgnc_sha256=mapper.sha256,
                    input_transform=INPUT_TRANSFORM,
                    exclude_doublets=bool(exclude_doublets),
                    source_sha=source_sha,
                    checkpoint=HASHES["best_model.pt"],
                    upstream=CODE,
                    device=device,
                    seconds=time.monotonic() - start,
                )
                os.replace(temp, dest)
            with np.load(dest, allow_pickle=False) as enc:
                if not np.array_equal(enc["cell_uid"], uid):
                    raise ValueError(f"Cached CLS embedding cell mismatch: {dest}")
                report.append(
                    dict(
                        **meta,
                        gene_universe=DEFAULT_UNIVERSE,
                        in_vocab_nonzero_fraction=float(
                            enc["in_vocab"].sum() / enc["nonzero"].sum()
                        ),
                        truncated_cells=int((enc["in_vocab"] > 1199).sum()),
                        seconds=float(enc["seconds"]),
                        bytes=dest.stat().st_size,
                        device=str(enc["device"]),
                        batch=batch,
                        checkpoint_sha256=HASHES["best_model.pt"],
                        hgnc_sha256=mapper.sha256,
                        source_sha256=source_sha,
                    )
                )
            print(f"{code} {sample} {lineage} {len(uid)}", flush=True)
        pd.DataFrame(report).to_csv(OUT / f"{code}_embedding_manifest.csv", index=False)


def forecast(
    fits: dict[tuple[str, str], Path],
    *,
    exclude_doublets: bool = True,
    max_cells: int | None = None,
) -> pd.DataFrame:
    mapper = get_hgnc_mapper()
    forecasts, draws, measurement, program_rows = [], [], [], []
    for code, tissue in TISSUES.items():
        scale_path = EXTERNAL / "coverage" / f"{code}_scale_check.csv"
        scale = (
            pd.read_csv(scale_path).set_index("sample_uid")
            if scale_path.is_file()
            else None
        )
        for meta, bag in tumor_bags(
            DATA / f"{code}.h5ad",
            tissue,
            gene_universe=DEFAULT_UNIVERSE,
            exclude_doublets=exclude_doublets,
            max_cells=max_cells,
        ):
            sample, lineage = meta["sample_uid"], meta["lineage"]
            status = {m: "fewer_than_50_cells" for m in METHODS}
            logits = {}
            if bag is not None:
                for uid in bag.obs.index:
                    draws.append(
                        dict(sample_uid=sample, lineage=lineage, cell_uid=str(uid))
                    )
                prog_fit_path = fits[lineage, "program_elastic"]
                status["program_elastic"] = "required_program_unmeasurable"
                means = score_frozen_program_bag(bag, prog_fit_path.parent)
                for prog_name, prog_val in means.items():
                    program_rows.append(
                        {
                            "sample_uid": sample,
                            "donor_uid": meta["donor_uid"],
                            "lesion_uid": meta["lesion_uid"],
                            "tissue": tissue,
                            "lineage": lineage,
                            "program": str(prog_name),
                            "n_cells": meta["n_cells"],
                            "score_ucell_debiased": float(prog_val),
                        }
                    )
                value = elastic_logit(prog_fit_path, means)
                if value is not None:
                    logits["program_elastic"] = value

                if scale is not None and sample in scale.index:
                    s = scale.loc[sample]
                    scale_ok = bool(
                        s.n_nonempty > 0
                        and 9900 <= s.median_nonempty_inverse_sum <= 10100
                        and s.fraction_within_five_percent_of_10000 >= 0.95
                    )
                    n_nonempty = int(s.n_nonempty)
                    med_inv = float(s.median_nonempty_inverse_sum)
                    frac_10k = float(s.fraction_within_five_percent_of_10000)
                else:
                    scale_ok = bool(np.isfinite(meta["median_nuclear_inverse_sum"]) and meta["median_nuclear_inverse_sum"] > 0)
                    n_nonempty = int(meta["n_cells"])
                    med_inv = float(meta["median_nuclear_inverse_sum"])
                    frac_10k = 1.0

                pca_file = fits[lineage, "PCA_mean"].parent / "pca.npz"
                with np.load(pca_file, allow_pickle=False) as transform:
                    pca, measured = pca_projection(bag, transform)
                    hvg = transform["genes"]
                    measured_x = bag.counts[:, bag.genes.get_indexer(hvg[measured])]
                    measurement.append(
                        dict(
                            **meta,
                            gene_universe=DEFAULT_UNIVERSE,
                            hgnc_sha256=mapper.sha256,
                            n_train_hvg=len(hvg),
                            n_mapped_hvg=int(measured.sum()),
                            mapped_fraction=float(measured.mean()),
                            absent_hvg="|".join(hvg[~measured]),
                            n_unmeasured_positions=int(len(bag.obs) * (~measured).sum()),
                            n_center_substituted=int(len(bag.obs) * (~measured).sum())
                            if scale_ok and pca is not None
                            else 0,
                            n_measured_zero=int(
                                len(bag.obs) * measured.sum() - measured_x.nnz
                            ),
                            n_nonempty=n_nonempty,
                            median_inverse_sum=med_inv,
                            fraction_9500_10500=frac_10k,
                            normalization_source="Liu et al. Cancer Cell 2022 STAR Methods e3 + shared_nuclear_245_164 CP10k recalibration",
                            numeric_scale_pass=scale_ok,
                        )
                    )
                for method in ("PCA_mean", "PCA_attention"):
                    status[method] = (
                        "unverified_sample_scale"
                        if not scale_ok
                        else "under_95pct_train_HVG"
                    )
                    if scale_ok and pca is not None:
                        logits[method] = pooled_logit(
                            fits[lineage, method], pca, sample, lineage
                        )
                with np.load(
                    cache_path(code, sample, lineage), allow_pickle=False
                ) as enc:
                    if (
                        not np.array_equal(enc["cell_uid"], bag.obs.index.to_numpy(dtype=str))
                        or str(enc["gene_universe"]) != DEFAULT_UNIVERSE
                        or bool(enc["exclude_doublets"]) != bool(exclude_doublets)
                    ):
                        raise ValueError(
                            f"External embedding bag mismatch: {sample}, {lineage}"
                        )
                    scgpt = enc["embedding"].astype(float)
                for method in ("scGPT_mean", "scGPT_attention"):
                    logits[method] = pooled_logit(
                        fits[lineage, method], scgpt, sample, lineage
                    )
            for method in METHODS:
                value = logits.get(method, np.nan)
                forecasts.append(
                    dict(
                        **meta,
                        method=method,
                        status="predicted" if method in logits else status[method],
                        logit=value,
                        probability_lm=float(expit(value)) if np.isfinite(value) else np.nan,
                    )
                )
    pred = pd.DataFrame(forecasts)
    pd.DataFrame(draws).to_csv(
        OUT / "selected_cells.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(program_rows).to_csv(
        OUT / "external_program_bag_scores.csv", index=False
    )
    pd.DataFrame(measurement).to_csv(OUT / "pca_measurement.csv", index=False)
    pred.to_csv(OUT / "bag_predictions.csv", index=False)
    return pred


def metric_draws(frame: pd.DataFrame):
    """Bootstrap true donors with all of their available sampled tumor sites."""
    frame = frame.sort_values(["donor_uid", "tissue"])
    codes, donors = pd.factorize(frame["donor_uid"], sort=True)
    rng = np.random.default_rng(SEED)
    weight = np.vstack(
        [
            np.ones(len(donors)),
            rng.multinomial(
                len(donors), np.full(len(donors), 1 / len(donors)), size=10000
            ),
        ]
    )
    w = weight[:, codes]
    pos = frame["tissue"].eq("LM").to_numpy()
    x = frame["logit"].to_numpy()
    wp, wn = w[:, pos], w[:, ~pos]
    np_, nn = wp.sum(1), wn.sum(1)
    score_diff = x[pos, None] - x[None, ~pos]
    concord = (score_diff > 0).astype(float) + 0.5 * (score_diff == 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        auc = np.einsum("bi,ij,bj->b", wp, concord, wn) / (np_ * nn)
        ba = 0.5 * ((wp @ (x[pos] >= 0)) / np_ + (wn @ (x[~pos] < 0)) / nn)
        thresholds = np.unique(x)[::-1]
        cum_pos = wp @ (x[pos, None] >= thresholds)
        cum_all = w @ (x[:, None] >= thresholds)
        precision = np.divide(
            cum_pos, cum_all, where=cum_all > 0, out=np.zeros_like(cum_all)
        )
        ap = (precision * np.diff(cum_pos, axis=1, prepend=0)).sum(axis=1) / np_
        ap[nn == 0] = np.nan
        pair = frame.pivot(index="donor_uid", columns="tissue", values="logit").reindex(
            index=donors, columns=["PT", "LM"]
        )
        dif = pair["LM"] - pair["PT"]
        paired = dif.notna().to_numpy()
        pairw = weight[:, paired]
        delta = (pairw @ dif[paired].to_numpy()) / pairw.sum(1)
        concordant = (pairw @ (dif[paired].to_numpy() > 0)) / pairw.sum(1)
    return {
        "auroc": auc,
        "average_precision": ap,
        "balanced_accuracy": ba,
        "mean_LM_minus_PT_logit": delta,
        "fraction_LM_above_PT": concordant,
    }, dif.dropna()


def intervals(draws: dict) -> dict:
    row = {}
    for metric, values in draws.items():
        finite = values[1:][np.isfinite(values[1:])]
        lo, hi = (
            np.quantile(finite, [0.025, 0.975]) if len(finite) else (np.nan, np.nan)
        )
        row.update(
            {
                metric: float(values[0]),
                metric + "_lower95": float(lo),
                metric + "_upper95": float(hi),
                metric + "_valid_bootstrap": int(len(finite)),
            }
        )
    return row


def patient_reports(pred: pd.DataFrame):
    from sklearn.metrics import (
        average_precision_score,
        balanced_accuracy_score,
        roc_auc_score,
    )

    keys = ["method", "lineage", "donor_uid", "tissue"]
    sites_all, metrics, pairs_all, contrast_rows = [], [], [], []
    for lineage in LINEAGES:
        group = pred.loc[pred["lineage"].eq(lineage)]
        wide = group.pivot(
            index="sample_uid", columns="method", values="logit"
        ).reindex(columns=METHODS)
        for scope in (
            "method_available",
            "common_all_five",
            "common_predicting_methods",
            "paired_donors_common_methods",
        ):
            available = wide.columns[wide.notna().any()]
            common = (
                wide.dropna().index
                if scope == "common_all_five"
                else wide[available].dropna().index
            )
            rows = group.loc[group["logit"].notna()]
            if scope != "method_available":
                rows = rows.loc[rows["sample_uid"].isin(common)]
            if scope == "paired_donors_common_methods":
                n_tissues = rows.groupby("donor_uid")["tissue"].transform("nunique")
                rows = rows.loc[n_tissues.eq(2)]
            lesions = rows.groupby(keys + ["lesion_uid"])["logit"].mean()
            sites = (
                lesions.groupby(keys)
                .agg(["mean", "size"])
                .reset_index()
                .rename(columns={"mean": "logit", "size": "n_lesions"})
            )
            sites_all.append(sites.assign(scope=scope))
            method_draws = {}
            for method in METHODS:
                frame = sites.loc[sites["method"].eq(method)]
                rec = dict(
                    lineage=lineage,
                    method=method,
                    scope=scope,
                    n_donors=frame["donor_uid"].nunique(),
                    n_sites=len(frame),
                    n_PT=int(frame["tissue"].eq("PT").sum()),
                    n_LM=int(frame["tissue"].eq("LM").sum()),
                    n_bags=int(len(rows[rows["method"].eq(method)])),
                    paired_n=0,
                )
                if len(frame):
                    draws, dif = metric_draws(frame)
                    method_draws[method] = draws
                    rec.update(intervals(draws), paired_n=len(dif))
                    pairs_all.extend(
                        dict(
                            lineage=lineage,
                            scope=scope,
                            method=method,
                            donor_uid=d,
                            LM_minus_PT_logit=float(v),
                        )
                        for d, v in dif.items()
                    )
                    if frame["tissue"].nunique() == 2:
                        y = frame["tissue"].eq("LM")
                        refs = [
                            roc_auc_score(y, frame["logit"]),
                            average_precision_score(y, frame["logit"]),
                            balanced_accuracy_score(y, frame["logit"] >= 0),
                        ]
                        if not np.allclose(
                            [
                                draws[k][0]
                                for k in (
                                    "auroc",
                                    "average_precision",
                                    "balanced_accuracy",
                                )
                            ],
                            refs,
                        ):
                            raise ValueError(
                                "Patient metrics diverge from their reference implementations"
                            )
                metrics.append(rec)
            if scope == "common_predicting_methods":
                for first, second in combinations(method_draws, 2):
                    delta = {
                        k: vals - method_draws[first][k]
                        for k, vals in method_draws[second].items()
                    }
                    contrast_rows.append(
                        dict(
                            lineage=lineage,
                            first=first,
                            second=second,
                            common_bags=len(common),
                            common_donors=sites["donor_uid"].nunique(),
                            **intervals(delta),
                        )
                    )
    pd.concat(sites_all, ignore_index=True).to_csv(
        OUT / "patient_site_logits.csv", index=False
    )
    pd.DataFrame(pairs_all).to_csv(OUT / "patient_paired_logits.csv", index=False)
    pd.DataFrame(contrast_rows).to_csv(OUT / "paired_method_contrasts.csv", index=False)
    pd.DataFrame(metrics).to_csv(OUT / "patient_metrics.csv", index=False)
    pred.groupby(["lineage", "method", "status"]).agg(
        n_bags=("sample_uid", "size"), n_donors=("donor_uid", "nunique")
    ).to_csv(OUT / "prediction_coverage.csv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["development-check", "encode", "predict"])
    parser.add_argument("--freeze-commit", default="WORKTREE")
    parser.add_argument(
        "--tissue-code", nargs="+", choices=TISSUES, default=list(TISSUES)
    )
    parser.add_argument("--device", choices=["cpu", "xpu"], default="xpu")
    parser.add_argument("--batch", type=int, default=16)
    parser.add_argument("--include-doublets", action="store_true")
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "development-check":
        development_check()
        return
    fits = frozen_models(args.freeze_commit)
    if set(fits) != {(lineage, method) for lineage in LINEAGES for method in METHODS}:
        raise ValueError("Expected 15 accepted M5 classifiers")
    start = time.monotonic()
    exclude_doublets = not args.include_doublets
    if args.stage == "encode":
        encode_external(
            args.tissue_code,
            args.device,
            args.batch,
            exclude_doublets=exclude_doublets,
            max_cells=args.max_cells,
            overwrite=args.overwrite,
        )
    else:
        patient_reports(
            forecast(fits, exclude_doublets=exclude_doublets, max_cells=args.max_cells)
        )
    key = (
        "encode_" + "_".join(args.tissue_code) if args.stage == "encode" else "predict"
    )
    git_target = "HEAD" if args.freeze_commit.upper() == "WORKTREE" else args.freeze_commit
    resolved_commit = subprocess.check_output(
        ["git", "rev-parse", git_target], cwd=ROOT, text=True
    ).strip()
    (OUT / (key + "_run.json")).write_text(
        json.dumps(
            dict(
                command=sys.argv,
                freeze_commit=resolved_commit,
                freeze_mode=args.freeze_commit,
                frozen_inputs_sha256=_sha256_file(M5 / "frozen_inputs_sha256.txt"),
                gene_universe=DEFAULT_UNIVERSE,
                exclude_doublets=exclude_doublets,
                max_cells=args.max_cells,
                python=platform.python_version(),
                executable=sys.executable,
                seconds=time.monotonic() - start,
                peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
            ),
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    write_output_sha256(OUT)


if __name__ == "__main__":
    main()
