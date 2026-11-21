"""Apply the frozen GSE245552 tissue classifiers to GSE164522 tumor specimens."""

import argparse
import hashlib
import json
import platform
from pathlib import Path
import resource
import subprocess
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]
from crlm.model_bags import Bag, LINEAGES, SEED, cell_seed, sample_bags
from linear_attention_mil import load_model, predict_bags


DATA = ROOT / "data/processed/m6/GSE164522"
EXTERNAL = ROOT / "artifacts/m6/GSE164522"
OUT = EXTERNAL / "model_predictions"
CACHE = ROOT / "embeddings/m6/GSE164522/scgpt"
M5 = ROOT / "artifacts/m5"
METHODS = (
    "program_elastic",
    "PCA_mean",
    "PCA_attention",
    "scGPT_mean",
    "scGPT_attention",
)
TISSUES = {"PT": "PT", "MT": "LM"}


def frozen_models(revision):
    frozen = M5 / "frozen_inputs_sha256.txt"
    version = subprocess.check_output(
        ["git", "show", f"{revision}:{frozen.relative_to(ROOT)}"], cwd=ROOT
    )
    if version != frozen.read_bytes():
        raise ValueError(
            "The freeze revision does not match the classifier input registry"
        )
    for line in version.decode().splitlines():
        expected, filename = line.split("  ", 1)
        if filename.startswith(
            (
                "models/m4/full/",
                "artifacts/m5/",
                "scripts/encode_scgpt.py",
                "src/crlm/model_bags.py",
            )
        ):
            with (ROOT / filename).open("rb") as source:
                if hashlib.file_digest(source, "sha256").hexdigest() != expected:
                    raise ValueError(f"Frozen input differs: {filename}")
    fits = pd.read_csv(M5 / "frozen_models.csv")
    return {
        (r.lineage, r.method): ROOT / r.path for r in fits.itertuples() if r.accepted
    }


def selected_rows(ids, rows, sample, lineage):
    if len(rows) < 50:
        return np.array([], dtype=int)
    if len(rows) <= 256:
        return rows
    ordered = rows[np.argsort(ids[rows].to_numpy())]
    rng = np.random.default_rng(cell_seed(f"{sample}|{lineage}"))
    return np.sort(rng.choice(ordered, 256, replace=False))


def tumor_bags(path, tissue):
    obj = ad.read_h5ad(path, backed="r")
    try:
        obs = obj.obs
        if set(obs.tissue) != {tissue} or not obs.index.is_unique:
            raise ValueError(f"Inconsistent cell identity or tissue in {path}")
        keep = obj.var.unambiguous_symbol.to_numpy(bool)
        genes = pd.Index(obj.var.loc[keep, "gene_symbol"])
        if not genes.is_unique:
            raise ValueError("Ambiguous gene symbol after filtering")
        for sample, idx in obs.groupby("sample_uid", observed=True).indices.items():
            sample_frame = obs.iloc[idx]
            if sample_frame[["donor_uid", "lesion_uid"]].nunique().max() != 1:
                raise ValueError(f"Unresolved specimen identity: {sample}")
            for lineage in LINEAGES:
                rows = idx[sample_frame.lineage.eq(lineage).to_numpy()]
                sel = selected_rows(obs.index, rows, sample, lineage)
                b = None
                if len(sel):
                    mat = obj.X[sel, :].tocsr()[:, keep]
                    mat.eliminate_zeros()
                    if np.any(mat.data < 0) or not np.isfinite(mat.data).all():
                        raise ValueError(
                            f"Expected nonnegative processed expression in {sample}"
                        )
                    b = Bag(obs.iloc[sel].copy(), mat, genes, len(rows))
                yield (
                    dict(
                        sample_uid=sample,
                        donor_uid=sample_frame.donor_uid.iloc[0],
                        lesion_uid=sample_frame.lesion_uid.iloc[0],
                        tissue=tissue,
                        lineage=lineage,
                        n_population=len(rows),
                        n_cells=len(sel),
                        n_captured=len(idx),
                    ),
                    b,
                )
    finally:
        obj.file.close()


def pca_projection(bag, transform):
    positions = bag.genes.get_indexer(transform["genes"])
    measured = positions >= 0
    if measured.mean() < 0.95:
        return None, measured
    # Missing features have their training means; centered contributions are exactly zero.
    values = bag.counts[:, positions[measured]].toarray().astype(np.float32)
    values -= transform["center"][measured]
    return (values @ transform["components"][:, measured].T).astype(float), measured


def elastic_logit(path, means):
    with np.load(path) as fit:
        values = means.reindex(fit["columns"]).to_numpy(float)
        if not np.isfinite(values).all():
            return None
        scaled = (values - fit["center"]) / fit["scale"]
        return float((scaled @ fit["coef"].T + fit["bias"])[0])


def pooled_logit(path, vector, sample, lineage):
    model, info = load_model(path)
    if info["lineage"] != lineage:
        raise ValueError(f"Classifier lineage differs from bag: {path}")
    prediction, _ = predict_bags(model, info, [vector], [sample])
    return float(prediction.logit.iloc[0])


def development_check():
    from crlm.program_scores import rank_scores

    design = pd.read_csv(ROOT / "artifacts/m4/oof_predictions.csv")
    base = design.loc[design.method.isin(METHODS) & design.seed.eq(SEED)]
    genes = pd.read_csv(
        ROOT / "artifacts/m3/program_definitions_draft/program_genes.csv"
    )
    genes = genes.loc[genes.on_discovery_axis & genes.role.eq("complete_program")]
    rows = []
    for lineage in LINEAGES:
        local = base.loc[base.lineage.eq(lineage)]
        ids = local.groupby("sample_uid").method.nunique()
        sample = min(ids.index[ids.eq(len(METHODS))])
        expected = local.loc[local.sample_uid.eq(sample)].set_index("method")
        fold = expected.fold.iloc[0]
        folder = ROOT / "models/m4" / lineage.replace("/", "_") / f"fold{fold}"
        path = ROOT / "data/processed/m2" / f"{sample}.h5ad"
        bag = sample_bags(path)[lineage]
        obj = ad.read_h5ad(path, backed="r")
        try:
            orig = np.flatnonzero(obj.obs.qc_retained & obj.obs.lineage_m2.eq(lineage))
            sel = selected_rows(obj.obs.index, orig, sample, lineage)
            if not obj.obs.index[sel].equals(bag.obs.index):
                raise ValueError(f"Different selected development cells in {sample}")
        finally:
            obj.file.close()
        defs = (
            genes[genes.lineage.eq(lineage)]
            .groupby("program")
            .gene_symbol.apply(list)
            .to_dict()
        )
        scores = rank_scores(bag.counts, bag.genes, defs).mean()
        mat = bag.counts.astype(np.float32)
        mat.data *= np.repeat(
            10000 / np.asarray(mat.sum(axis=1)).ravel(), np.diff(mat.indptr)
        )
        np.log1p(mat.data, out=mat.data)
        with np.load(folder / "pca.npz") as transform:
            pca, _ = pca_projection(
                Bag(bag.obs, mat, bag.genes, bag.population), transform
            )
        with np.load(ROOT / f"embeddings/m4/scgpt/{sample}.npz") as encoded:
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


def cache_path(code, sample, lineage):
    return CACHE / code / (sample + "_" + lineage.replace("/", "_") + ".npz")


def encode_external(codes, device, batch):
    import torch
    from encode_scgpt import CODE, HASHES, encode, load_encoder, tokens

    profile = json.loads(
        (ROOT / "artifacts/m4/resource_profiles/scgpt_profile.json").read_text()
    )
    if profile["checkpoint"] != HASHES or not profile["devices"][device]["accepted"]:
        raise ValueError(f"Unprofiled checkpoint/device: {device}")
    torch.set_num_threads(4)
    model, vocab = load_encoder()
    model.to(device)
    for code in codes:
        path = DATA / f"{code}.h5ad"
        (CACHE / code).mkdir(parents=True, exist_ok=True)
        obj = ad.read_h5ad(path, backed="r")
        source_sha = obj.uns["expression_sha256"]
        obj.file.close()
        report = []
        for meta, bag in tumor_bags(path, TISSUES[code]):
            if bag is None:
                continue
            sample, lineage = meta["sample_uid"], meta["lineage"]
            dest = cache_path(code, sample, lineage)
            uid = bag.obs.index.to_numpy(str)
            if not dest.exists():
                start = time.monotonic()
                cells = tokens(bag, vocab)
                embedded = encode(model, cells, vocab, batch, device)
                temp = dest.with_suffix(".partial.npz")
                np.savez_compressed(
                    temp,
                    embedding=embedded,
                    cell_uid=uid,
                    nonzero=np.array([x[3] for x in cells]),
                    in_vocab=np.array([x[4] for x in cells]),
                    tokens=np.array([len(x[1]) for x in cells]),
                    source_sha=source_sha,
                    checkpoint=HASHES["best_model.pt"],
                    upstream=CODE,
                    device=device,
                    seconds=time.monotonic() - start,
                )
                temp.replace(dest)
            with np.load(dest) as enc:
                if (
                    not np.array_equal(enc["cell_uid"], uid)
                    or enc["source_sha"].item() != source_sha
                    or enc["checkpoint"].item() != HASHES["best_model.pt"]
                    or enc["upstream"].item() != CODE
                ):
                    raise ValueError(f"Cached CLS embedding identity differs: {dest}")
                report.append(
                    dict(
                        **meta,
                        in_vocab_nonzero_fraction=enc["in_vocab"].sum()
                        / enc["nonzero"].sum(),
                        truncated_cells=int((enc["in_vocab"] > 1199).sum()),
                        seconds=enc["seconds"].item(),
                        bytes=dest.stat().st_size,
                        device=enc["device"].item(),
                        batch=batch,
                        checkpoint_sha256=HASHES["best_model.pt"],
                        source_sha256=source_sha,
                    )
                )
            print(f"{code} {sample} {lineage} {len(uid)}", flush=True)
        pd.DataFrame(report).to_csv(OUT / f"{code}_embedding_manifest.csv", index=False)


def forecast(fits):
    inputs = pd.read_csv(OUT / "external_model_inputs.csv")
    program_bags = dict(tuple(inputs.groupby(["sample_uid", "lineage"])))
    scores = pd.read_csv(
        DATA / "frozen_cell_scores.csv.gz",
        usecols=["sample_uid", "lineage", "program", "cell_uid", "score"],
    )
    scores = scores.loc[scores.sample_uid.isin(inputs.sample_uid)]
    programs = pd.read_csv(M5 / "gene_sets/model_program_genes.csv")
    cols = programs.groupby("lineage").program.unique()
    forecasts, draws, measurement = [], [], []
    for code, tissue in TISSUES.items():
        scale = pd.read_csv(
            EXTERNAL / "coverage" / f"{code}_scale_check.csv"
        ).set_index("sample_uid")
        for meta, bag in tumor_bags(DATA / f"{code}.h5ad", tissue):
            sample, lineage = meta["sample_uid"], meta["lineage"]
            status = {m: "fewer_than_50_cells" for m in METHODS}
            logits = {}
            if bag is not None:
                for uid in bag.obs.index:
                    draws.append(dict(sample_uid=sample, lineage=lineage, cell_uid=uid))
                program_data = program_bags.get((sample, lineage))
                status["program_elastic"] = "required_program_unmeasurable"
                if program_data is not None:
                    summary = program_data.set_index("program")
                    if (
                        summary.index.has_duplicates
                        or not (summary.n_cells == meta["n_cells"]).all()
                    ):
                        raise ValueError(
                            f"Program sample bag mismatch: {sample}, {lineage}"
                        )
                    cell = scores.loc[
                        scores.sample_uid.eq(sample)
                        & scores.lineage.eq(lineage)
                        & scores.program.isin(cols[lineage])
                        & scores.cell_uid.isin(bag.obs.index)
                    ]
                    counts = (
                        cell.groupby("program")
                        .cell_uid.nunique()
                        .reindex(summary.index)
                    )
                    means = cell.groupby("program").score.mean().reindex(summary.index)
                    if not counts.eq(meta["n_cells"]).all() or not np.allclose(
                        means, summary.score, atol=1e-12
                    ):
                        raise ValueError(
                            f"Frozen program scores use different selected cells: {sample}, {lineage}"
                        )
                    value = elastic_logit(fits[lineage, "program_elastic"], means)
                    if value is not None:
                        logits["program_elastic"] = value
                s = scale.loc[sample]
                scale_ok = (
                    s.n_nonempty > 0
                    and 9900 <= s.median_nonempty_inverse_sum <= 10100
                    and s.fraction_within_five_percent_of_10000 >= 0.95
                )
                pca_file = fits[lineage, "PCA_mean"].parent / "pca.npz"
                with np.load(pca_file) as transform:
                    pca, measured = pca_projection(bag, transform)
                    hvg = transform["genes"]
                    measured_x = bag.counts[:, bag.genes.get_indexer(hvg[measured])]
                    measurement.append(
                        dict(
                            **meta,
                            n_train_hvg=len(hvg),
                            n_mapped_hvg=int(measured.sum()),
                            mapped_fraction=measured.mean(),
                            absent_hvg="|".join(hvg[~measured]),
                            n_unmeasured_positions=int(
                                len(bag.obs) * (~measured).sum()
                            ),
                            n_center_substituted=int(len(bag.obs) * (~measured).sum())
                            if scale_ok and pca is not None
                            else 0,
                            n_measured_zero=int(
                                len(bag.obs) * measured.sum() - measured_x.nnz
                            ),
                            n_nonempty=s.n_nonempty,
                            median_inverse_sum=s.median_nonempty_inverse_sum,
                            fraction_9500_10500=s.fraction_within_five_percent_of_10000,
                            normalization_source="Liu et al. Cancer Cell 2022, STAR Methods e3",
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
                with np.load(cache_path(code, sample, lineage)) as enc:
                    if not np.array_equal(enc["cell_uid"], bag.obs.index.to_numpy(str)):
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
                        probability_lm=expit(value),
                    )
                )
    pred = pd.DataFrame(forecasts)
    pd.DataFrame(draws).to_csv(
        OUT / "selected_cells.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame(measurement).to_csv(OUT / "pca_measurement.csv", index=False)
    pred.to_csv(OUT / "bag_predictions.csv", index=False)
    return pred


def metric_draws(frame):
    """Bootstrap true donors with all of their available sampled tumor sites."""
    frame = frame.sort_values(["donor_uid", "tissue"])
    codes, donors = pd.factorize(frame.donor_uid, sort=True)
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
    pos = frame.tissue.eq("LM").to_numpy()
    x = frame.logit.to_numpy()
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
        dif = pair.LM - pair.PT
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


def intervals(draws):
    row = {}
    for metric, values in draws.items():
        finite = values[1:][np.isfinite(values[1:])]
        lo, hi = (
            np.quantile(finite, [0.025, 0.975]) if len(finite) else (np.nan, np.nan)
        )
        row.update(
            {
                metric: values[0],
                metric + "_lower95": lo,
                metric + "_upper95": hi,
                metric + "_valid_bootstrap": len(finite),
            }
        )
    return row


def patient_reports(pred):
    from itertools import combinations
    from sklearn.metrics import (
        average_precision_score,
        balanced_accuracy_score,
        roc_auc_score,
    )

    keys = ["method", "lineage", "donor_uid", "tissue"]
    sites_all, metrics, pairs_all, contrast_rows = [], [], [], []
    for lineage in LINEAGES:
        group = pred.loc[pred.lineage.eq(lineage)]
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
            rows = group.loc[group.logit.notna()]
            if scope != "method_available":
                rows = rows.loc[rows.sample_uid.isin(common)]
            if scope == "paired_donors_common_methods":
                n_tissues = rows.groupby("donor_uid").tissue.transform("nunique")
                rows = rows.loc[n_tissues.eq(2)]
            lesions = rows.groupby(keys + ["lesion_uid"]).logit.mean()
            sites = (
                lesions.groupby(keys)
                .agg(["mean", "size"])
                .reset_index()
                .rename(columns={"mean": "logit", "size": "n_lesions"})
            )
            sites_all.append(sites.assign(scope=scope))
            method_draws = {}
            for method in METHODS:
                frame = sites.loc[sites.method.eq(method)]
                rec = dict(
                    lineage=lineage,
                    method=method,
                    scope=scope,
                    n_donors=frame.donor_uid.nunique(),
                    n_sites=len(frame),
                    n_PT=int(frame.tissue.eq("PT").sum()),
                    n_LM=int(frame.tissue.eq("LM").sum()),
                    n_bags=int(len(rows[rows.method.eq(method)])),
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
                            LM_minus_PT_logit=v,
                        )
                        for d, v in dif.items()
                    )
                    if frame.tissue.nunique() == 2:
                        y = frame.tissue.eq("LM")
                        refs = [
                            roc_auc_score(y, frame.logit),
                            average_precision_score(y, frame.logit),
                            balanced_accuracy_score(y, frame.logit >= 0),
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
                            common_donors=sites.donor_uid.nunique(),
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
    parser.add_argument("--freeze-commit")
    parser.add_argument(
        "--tissue-code", nargs="+", choices=TISSUES, default=list(TISSUES)
    )
    parser.add_argument("--device", choices=["cpu", "xpu"], default="xpu")
    parser.add_argument("--batch", type=int, default=4)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "development-check":
        development_check()
        return
    if not args.freeze_commit:
        parser.error("The frozen M5 revision is required for external classification")
    fits = frozen_models(args.freeze_commit)
    if set(fits) != {(lineage, method) for lineage in LINEAGES for method in METHODS}:
        raise ValueError("Expected 15 accepted M5 classifiers")
    start = time.monotonic()
    if args.stage == "encode":
        encode_external(args.tissue_code, args.device, args.batch)
    else:
        patient_reports(forecast(fits))
    key = (
        "encode_" + "_".join(args.tissue_code) if args.stage == "encode" else "predict"
    )
    (OUT / (key + "_run.json")).write_text(
        json.dumps(
            dict(
                command=sys.argv,
                freeze_commit=subprocess.check_output(
                    ["git", "rev-parse", args.freeze_commit], cwd=ROOT, text=True
                ).strip(),
                python=platform.python_version(),
                executable=sys.executable,
                seconds=time.monotonic() - start,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                * 1024,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
