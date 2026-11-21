"""Combine OOF runs, verify patient splits and quantify held-out attribution stability."""

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

from benchmark_discovery import (
    ROOT,
    OUT,
    SPLIT,
    SEED,
    MODELS,
    library_inputs,
    scgpt_rep,
    bag_forward,
    unpack,
    summarize,
)
from linear_attention_mil import load_model


def scores(model, x):
    z = (x - model["center"]) / model["scale"]
    w, b, v = unpack(model["theta"], x.shape[1], "attention")
    _, _, _, a = bag_forward(z, w, b, v)
    return a, a * (z @ w)


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
    folds = json.loads(SPLIT.read_text())
    if (
        oof.duplicated(["method", "seed", "lineage", "sample_uid"]).any()
        or not np.array_equal(oof.fold, oof.donor_uid.map(folds["donor_test_fold"]))
        or not np.isfinite(oof.logit).all()
    ):
        raise ValueError("Invalid or duplicate OOF prediction")
    rows, groups = library_inputs(write_selection=False)
    split_rows = []
    stability = []
    for lineage, group in groups.items():
        sc = scgpt_rep(group)
        for fold in folds["folds"]:
            f = fold["fold"]
            dest = MODELS / lineage.replace("/", "_") / f"fold{f}"
            te = rows[
                rows.sample_uid.isin(group) & rows.donor_uid.isin(fold["test_donors"])
            ]
            for r in rows[rows.sample_uid.isin(group)].itertuples():
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
                for r in te.itertuples():
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
    contrasts.to_csv(OUT / "paired_model_comparisons/auc_differences.csv", index=False)
    manifest = []
    for path in sorted((ROOT / "embeddings/m4/scgpt").glob("GSM*.npz")):
        with np.load(path, allow_pickle=False) as d:
            h = d["embedding"]
            assert h.shape == (len(d["cell_uid"]), 512) and np.isfinite(h).all()
            assert np.allclose(np.linalg.norm(h, axis=1), 1, atol=2e-5)
            manifest.append(
                dict(
                    sample_uid=path.stem,
                    cells=len(h),
                    width=512,
                    seconds=float(d["seconds"]),
                    bytes=path.stat().st_size,
                    total_nonzero=int(d["expressed_counts"].sum()),
                    matched_nonzero=int(d["in_vocab"].sum()),
                    truncated_cells=int(np.sum(d["in_vocab"] > 1199)),
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                )
            )
    enc = pd.DataFrame(manifest)
    enc.to_csv(OUT / "resource_profiles/scgpt_cache_manifest.csv", index=False)
    selected = pd.read_csv(
        OUT / "selected_cells.csv.gz", usecols=["cell_uid", "sample_uid"]
    )
    assert selected.cell_uid.is_unique and len(selected) == enc.cells.sum()
    print(
        f"{len(runs)} fits, {len(oof)} OOF scores, {len(selected)} selected cells; elapsed {time.perf_counter() - start:.1f}s"
    )


if __name__ == "__main__":
    with threadpool_limits(limits=4):
        main()
