"""Refit the 15 prespecified transferable classifiers on full development tumor bags (`245_164_shared_nuclear`)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter"), str(ROOT / "scripts")]

from benchmark_discovery import (
    DEFAULT_UNIVERSE,
    MODELS,
    OUT,
    SEED,
    build_design,
    elastic,
    fit_model,
    library_inputs,
    load_lineage_programs,
    pca_rep,
    save_model,
    scgpt_rep,
)
from crlm.hgnc_alias_map import get_hgnc_mapper, load_shared_nuclear_universe


def _sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    mapper = get_hgnc_mapper()
    universe_genes = list(load_shared_nuclear_universe(DEFAULT_UNIVERSE))
    all_program_defs = load_lineage_programs(universe_genes)
    rows, groups = library_inputs(write_selection=False)
    models = []
    for lineage, group in groups.items():
        tr = rows[rows["sample_uid"].isin(group)].reset_index(drop=True)
        design = build_design(tr)
        dest = MODELS / "full" / lineage.replace("/", "_")
        dest.mkdir(parents=True, exist_ok=True)
        pca, _ = pca_rep(group, design.samples, dest / "pca.npz")
        sc = scgpt_rep(group)
        rank_path = dest / "program_elastic.npz"
        _, ok, iterations, nparams = elastic(
            group,
            tr,
            [],
            "program_elastic",
            rank_path,
            lineage=lineage,
            program_defs=all_program_defs[lineage],
            universe_genes=universe_genes,
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
                            training_donors=sorted(tr["donor_uid"].unique()),
                            representation=rep,
                            gene_universe=DEFAULT_UNIVERSE,
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
        n_cells_total = int(sum(len(group[s][0].obs) for s in design.samples))
        (dest / "training_identity.json").write_text(
            json.dumps(
                {
                    "lineage": lineage,
                    "donors": sorted(tr["donor_uid"].unique()),
                    "samples": design.samples,
                    "n_cells_total": n_cells_total,
                    "gene_universe": DEFAULT_UNIVERSE,
                    "n_universe_genes": len(universe_genes),
                    "hgnc_sha256": mapper.sha256,
                    "seed": SEED,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    table = pd.DataFrame(models)
    table["sha256"] = [
        _sha256_file(ROOT / p) if p else ""
        for p in table["path"]
    ]
    OUT.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT / "full_data_models.csv", index=False)

    # Refresh output_sha256.txt if present so full_data_models.csv is included
    sha_manifest = OUT / "output_sha256.txt"
    if sha_manifest.is_file():
        files = sorted(
            p
            for p in OUT.rglob("*")
            if p.is_file() and p.name not in ("output_sha256.txt", ".gitkeep")
        )
        sha_manifest.write_text(
            "".join(f"{_sha256_file(p)}  {p.relative_to(ROOT)}\n" for p in files),
            encoding="utf-8",
        )
    print(table.groupby(["method", "accepted"]).size().to_string())


if __name__ == "__main__":
    with threadpool_limits(limits=4):
        main()
