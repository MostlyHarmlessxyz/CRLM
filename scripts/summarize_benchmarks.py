"""Combine OOF runs, verify patient splits, summarize direct-prior & sub-bag benchmarks, and lock SHA256s."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter"), str(ROOT / "scripts")]

from benchmark_discovery import (
    DEFAULT_UNIVERSE,
    MODELS,
    OUT,
    SEED,
    SPLIT,
    bag_forward,
    library_inputs,
    scgpt_rep,
    summarize,
    unpack,
)
from crlm.hgnc_alias_map import load_shared_nuclear_universe
from crlm.model_bags import JOINT_LINEAGE, sample_joint_bag
from crlm.prior_mil_metrics import run as run_prior_mil_metrics
from linear_attention_mil import load_model


def _sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def scores(model, x):
    z = (x - model["center"]) / model["scale"]
    w, b, v = unpack(model["theta"], x.shape[1], "attention")
    _, _, _, a = bag_forward(z, w, b, v)
    return a, a * (z @ w)


def write_output_sha256(stage_dir: Path) -> int:
    target = stage_dir / "output_sha256.txt"
    files = sorted(
        p
        for p in stage_dir.rglob("*")
        if p.is_file() and p.name not in ("output_sha256.txt", ".gitkeep")
    )
    lines = [f"{_sha256_file(p)}  {p.relative_to(ROOT)}\n" for p in files]
    target.write_text("".join(lines), encoding="utf-8")
    return len(files)


def main():
    start = time.perf_counter()
    sources = [
        OUT / "PCA",
        OUT / "scGPT",
        OUT / "direct_prior/PCA32",
        OUT / "direct_prior/scGPT32",
        OUT / "direct_prior/scGPT512",
    ]
    oof = pd.concat(
        [pd.read_csv(p / "oof_predictions.csv") for p in sources], ignore_index=True
    )
    runs = pd.concat(
        [pd.read_csv(p / "benchmark_runs.csv") for p in sources], ignore_index=True
    )
    folds = json.loads(SPLIT.read_text(encoding="utf-8"))
    if (
        oof.duplicated(["method", "seed", "lineage", "sample_uid"]).any()
        or not np.array_equal(oof["fold"], oof["donor_uid"].map(folds["donor_test_fold"]))
        or not np.isfinite(oof["logit"]).all()
    ):
        raise ValueError("Invalid or duplicate OOF prediction")

    rows, groups = library_inputs(write_selection=False)
    universe_genes = list(load_shared_nuclear_universe(DEFAULT_UNIVERSE))
    joint_samples = set()
    joint_cell_uids: set[str] = set()
    for gsm in rows["sample_uid"]:
        jb = sample_joint_bag(
            ROOT / f"data/processed/m2/{gsm}.h5ad",
            min_cells=50,
            max_cells=None,
            gene_universe=universe_genes,
        )
        if jb is not None:
            joint_samples.add(gsm)
            joint_cell_uids.update(jb.obs.index.astype(str).tolist())

    split_rows = []
    stability = []
    for lineage, group in groups.items():
        sc = scgpt_rep(group)
        for fold in folds["folds"]:
            f = fold["fold"]
            dest = MODELS / lineage.replace("/", "_") / f"fold{f}"
            te = rows[
                rows["sample_uid"].isin(group) & rows["donor_uid"].isin(fold["test_donors"])
            ]
            for r in rows[rows["sample_uid"].isin(group)].itertuples(index=False):
                split_rows.append(
                    dict(
                        fold=f,
                        lineage=lineage,
                        donor_uid=r.donor_uid,
                        sample_uid=r.sample_uid,
                        tissue=r.tissue,
                        role="test" if r.donor_uid in fold["test_donors"] else "train",
                    )
                )
            with np.load(dest / "pca.npz", allow_pickle=False) as pp:
                hv, center, comp = pp["genes"], pp["center"], pp["components"]
            for rep in ("PCA", "scGPT"):
                fits = {}
                for seed in (SEED, SEED + 1, SEED + 2):
                    path = dest / f"{rep}_attention_{seed}.npz"
                    if path.exists():
                        fits[seed] = load_model(path)[0]
                if SEED not in fits:
                    continue
                for r in te.itertuples(index=False):
                    s = r.sample_uid
                    bag, x, _ = group[s]
                    if rep == "PCA":
                        ix = bag.genes.get_indexer(hv)
                        if min(ix) < 0:
                            raise ValueError("Missing development PCA feature")
                        x = (x[:, ix].toarray() - center) @ comp.T
                    else:
                        x = sc[s]
                    ref = scores(fits[SEED], x)
                    k = max(1, len(x) // 10)
                    for seed, model in fits.items():
                        if seed == SEED:
                            continue
                        other = scores(model, x)
                        top0 = set(np.argsort(-np.abs(ref[1]))[:k])
                        top1 = set(np.argsort(-np.abs(other[1]))[:k])
                        stability.append(
                            dict(
                                lineage=lineage,
                                rep=rep,
                                fold=f,
                                sample_uid=s,
                                other_seed=seed,
                                signed_rho=spearmanr(ref[1], other[1]).statistic,
                                attention_rho=spearmanr(ref[0], other[0]).statistic,
                                top_abs_jaccard=len(top0 & top1) / len(top0 | top1),
                            )
                        )

    for fold in folds["folds"]:
        f = fold["fold"]
        for r in rows[rows["sample_uid"].isin(joint_samples)].itertuples(index=False):
            split_rows.append(
                dict(
                    fold=f,
                    lineage=JOINT_LINEAGE,
                    donor_uid=r.donor_uid,
                    sample_uid=r.sample_uid,
                    tissue=r.tissue,
                    role="test" if r.donor_uid in fold["test_donors"] else "train",
                )
            )

    splits = OUT / "split_manifest"
    splits.mkdir(exist_ok=True)
    pd.DataFrame(split_rows).to_csv(
        splits / "eligible_samples_by_fold.csv", index=False
    )
    pd.DataFrame(stability).to_csv(OUT / "attribution_stability.csv", index=False)
    oof.to_csv(OUT / "oof_predictions.csv", index=False)
    runs.to_csv(OUT / "benchmark_runs.csv", index=False)

    patient, metrics, contrasts = summarize(oof, SEED)
    patient.to_csv(OUT / "patient_site_logits.csv", index=False)
    metrics.to_csv(OUT / "patient_metrics.csv", index=False)
    compare = OUT / "paired_model_comparisons"
    compare.mkdir(exist_ok=True)
    contrasts.to_csv(compare / "auc_differences.csv", index=False)

    # Generate direct_prior matched-mean comparisons across single lineages and All3Lineages
    run_prior_mil_metrics()

    # Combine sub-bag 256x10 variance sensitivity tables
    subbag_parts = [
        pd.read_csv(OUT / rep / "subbag_256x10_sensitivity.csv")
        for rep in ("PCA", "scGPT")
        if (OUT / rep / "subbag_256x10_sensitivity.csv").is_file()
    ]
    if subbag_parts:
        pd.concat(subbag_parts, ignore_index=True).to_csv(
            OUT / "subbag_256x10_sensitivity.csv", index=False
        )

    manifest = []
    encoded_uids: set[str] = set()
    for path in sorted((ROOT / "embeddings/m4/scgpt").glob("GSM*.npz")):
        with np.load(path, allow_pickle=False) as d:
            h = d["embedding"]
            uids = [str(u) for u in d["cell_uid"].tolist()]
            assert h.shape == (len(uids), 512) and np.isfinite(h).all()
            assert np.allclose(np.linalg.norm(h, axis=1), 1.0, atol=2e-5)
            encoded_uids.update(uids)
            manifest.append(
                dict(
                    sample_uid=path.stem,
                    cells=len(h),
                    width=512,
                    gene_universe=str(d["gene_universe"]) if "gene_universe" in d.files else "",
                    seconds=float(d["seconds"]),
                    bytes=path.stat().st_size,
                    total_nonzero=int(d["expressed_counts"].sum()),
                    matched_nonzero=int(d["in_vocab"].sum()),
                    truncated_cells=int(np.sum(d["in_vocab"] > 1199)),
                    sha256=_sha256_file(path),
                )
            )
    enc = pd.DataFrame(manifest)
    (OUT / "resource_profiles").mkdir(exist_ok=True)
    enc.to_csv(OUT / "resource_profiles/scgpt_cache_manifest.csv", index=False)

    selected = pd.read_csv(
        OUT / "selected_cells.csv.gz", usecols=["cell_uid", "sample_uid"]
    )
    sel_uids = set(selected["cell_uid"].astype(str).tolist())
    assert selected["cell_uid"].is_unique
    assert sel_uids.issubset(encoded_uids), "Selected single-lineage cells missing from scGPT cache"
    assert joint_cell_uids == encoded_uids, (
        f"Encoded scGPT cell set ({len(encoded_uids)}) differs from All3Lineages QC-retained set ({len(joint_cell_uids)})"
    )

    n_hashed = write_output_sha256(OUT)
    print(
        f"{len(runs)} fits, {len(oof)} OOF scores, {len(selected)} single-lineage cells, "
        f"{len(encoded_uids)} total encoded tumor cells, {n_hashed} files hashed; "
        f"elapsed {time.perf_counter() - start:.1f}s"
    )


if __name__ == "__main__":
    with threadpool_limits(limits=4):
        main()
