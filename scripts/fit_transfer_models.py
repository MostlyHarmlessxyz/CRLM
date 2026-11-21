"""Refit the prespecified transferable classifiers on eligible development bags."""

import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
from threadpoolctl import threadpool_limits

from benchmark_discovery import (
    ROOT,
    OUT,
    MODELS,
    SEED,
    build_design,
    elastic,
    fit_model,
)
from benchmark_discovery import library_inputs, pca_rep, save_model, scgpt_rep


def main():
    rows, groups = library_inputs(write_selection=False)
    models = []
    for lineage, group in groups.items():
        tr = rows[rows.sample_uid.isin(group)].reset_index(drop=True)
        design = build_design(tr)
        dest = MODELS / "full" / lineage.replace("/", "_")
        dest.mkdir(parents=True, exist_ok=True)
        pca = pca_rep(group, design.samples, dest / "pca.npz")
        sc = scgpt_rep(group)
        rank_path = dest / "program_elastic.npz"
        _, ok, iterations, nparams = elastic(
            group, tr, [], "program_elastic", rank_path
        )
        models.append(
            dict(
                lineage=lineage,
                method="program_elastic",
                accepted=ok,
                iterations=iterations,
                parameters=nparams,
                path=str(rank_path.relative_to(ROOT)),
            )
        )
        for rep, h in [("PCA", pca), ("scGPT", sc)]:
            for pooling in ("mean", "attention"):
                model, fit = fit_model(
                    [h[s] for s in design.samples],
                    design,
                    pooling=pooling,
                    lambda_pair=0.1,
                    lambda_l2=0.01,
                    seed=SEED,
                    standardize=True,
                    maxiter=500,
                )
                path = dest / f"{rep}_{pooling}.npz"
                if fit["converged"]:
                    save_model(
                        path,
                        model,
                        dict(
                            lineage=lineage,
                            fit=fit,
                            training_samples=design.samples,
                            training_donors=sorted(tr.donor_uid.unique()),
                            representation=rep,
                        ),
                    )
                models.append(
                    dict(
                        lineage=lineage,
                        method=f"{rep}_{pooling}",
                        accepted=fit["converged"],
                        parameters=fit["n_parameters"],
                        iterations=fit["iterations"],
                        path=str(path.relative_to(ROOT)) if fit["converged"] else "",
                    )
                )
        (dest / "training_identity.json").write_text(
            json.dumps(
                {
                    "lineage": lineage,
                    "donors": sorted(tr.donor_uid.unique()),
                    "samples": design.samples,
                    "seed": SEED,
                },
                indent=2,
            )
            + "\n"
        )
    table = pd.DataFrame(models)
    table["sha256"] = [
        hashlib.sha256((ROOT / p).read_bytes()).hexdigest() if p else ""
        for p in table.path
    ]
    table.to_csv(OUT / "full_data_models.csv", index=False)
    print(table.groupby(["method", "accepted"]).size().to_string())


if __name__ == "__main__":
    with threadpool_limits(limits=4):
        main()
