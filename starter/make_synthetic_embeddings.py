#!/usr/bin/env python3
"""Make tiny, explicitly synthetic cell embeddings for the MIL smoke demo."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from _common import read_csv, stable_seed, write_json


def generate(manifest, outdir, seed=20261121):
    frame = read_csv(manifest)
    bags = {}
    for row in frame.itertuples(index=False):
        rng = np.random.default_rng(stable_seed(seed, "bag", row.sample_uid))
        donor_rng = np.random.default_rng(stable_seed(seed, "donor", row.donor_uid))
        n_cells = int(rng.integers(17, 34))
        h = rng.normal(0, 0.7, (n_cells, 5)) + donor_rng.normal(0, 0.2, 5)
        # Deliberately manufactured tissue-associated signals, not real data.
        h[:, 0] += {"PT": -0.8, "LM": 0.8, "NC": -0.2, "NL": 0.5}[row.tissue]
        h[:, 1] += 0.4 if row.tissue in {"LM", "NL"} else -0.4
        if row.tissue == "LM":
            h[:max(1, n_cells // 4), 2] += 1.0
        bags[row.sample_uid] = h
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(outdir / "cell_embeddings.npz", **bags)
    provenance = {
        "lineage": "SyntheticImmune", "encoder_id": "synthetic_generator_not_scGPT_or_PCA",
        "representation_frozen": True, "preprocessing_fit_scope": "synthetic",
        "gene_universe_policy": "No genes: five manufactured numerical features. Never interpret as biological expression.",
        "gene_coverage_audited": True, "embedding_dimension": 5, "seed": int(seed),
        "note": "Uses tissue labels only to construct a toy example. Performance has no biological meaning.",
    }
    write_json(outdir / "embedding_provenance.json", provenance)
    return bags, provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(Path(__file__).parent / "synthetic_demo/sample_manifest.csv"))
    parser.add_argument("--outdir", default=str(Path(__file__).parent / "synthetic_demo"))
    parser.add_argument("--seed", type=int, default=20261121)
    args = parser.parse_args()
    bags, _ = generate(args.manifest, args.outdir, args.seed)
    print(f"Wrote {len(bags)} synthetic bags, five features each; this is not a biological result.")


if __name__ == "__main__":
    main()
