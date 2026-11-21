"""Select the same development cells for every representation and pooling method."""

from dataclasses import dataclass
from pathlib import Path
import hashlib

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261121
LINEAGES = ("T/NK", "Mono/Mac/DC", "B/Plasma")


def cell_seed(uid, seed=SEED):
    """Stable random stream for a cell or bag, independent of processing order."""
    return int.from_bytes(
        hashlib.sha256(f"{seed}|{uid}".encode()).digest()[:4], "little"
    )


@dataclass
class Bag:
    obs: pd.DataFrame
    counts: object
    genes: pd.Index
    population: int


def sample_bags(path):
    obj = ad.read_h5ad(path, backed="r")
    try:
        obs = obj.obs
        if obs.tissue.nunique() != 1 or obs.tissue.iloc[0] not in ("LM", "PT"):
            raise ValueError(f"Expected a tumor library at {path}")
        retained = obs.qc_retained.to_numpy(dtype=bool)
        if not obs.index.is_unique or not obj.var.gene_symbol.is_unique:
            raise ValueError(f"Non-unique cell or gene identifiers in {path}")
        groups = {}
        for lineage in LINEAGES:
            rows = np.flatnonzero(retained & obs.lineage_m2.eq(lineage).to_numpy())
            population = len(rows)
            if population < 50:
                continue
            if population > 256:
                ordered = rows[np.argsort(obs.index[rows].to_numpy())]
                rng = np.random.default_rng(cell_seed(f"{path.stem}|{lineage}"))
                rows = np.sort(rng.choice(ordered, 256, replace=False))
            frame = obs.iloc[rows][
                [
                    "donor_uid",
                    "lesion_uid",
                    "sample_uid",
                    "tissue",
                    "lineage_m2",
                    "subtype_m2",
                    "n_counts_visible",
                    "n_genes_visible",
                ]
            ].copy()
            groups[lineage] = Bag(
                frame, obj.X[rows, :].tocsr(), pd.Index(obj.var.gene_symbol), population
            )
        return groups
    finally:
        obj.file.close()
