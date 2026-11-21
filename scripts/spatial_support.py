"""Measure within-section CD3/CD8 proximity to SPP1-positive and negative myeloid cells."""

import gzip
import hashlib
from pathlib import Path
import re
import time

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/GSE335552"
OUT = ROOT / "artifacts/m7/spatial_measurement_audit"
PANEL = "CD68 CD163 C1QA CD14 LYZ EPCAM CD3D CD3E CD8A CD8B SPP1 CXCL9 CXCL10".split()


def spatial_difference(coords, distance, positive, tile_size, rng):
    tile = np.floor(coords / tile_size).astype(int)
    codes = np.unique(tile, axis=0, return_inverse=True)[1]
    chunks = []
    order = np.argsort(codes)
    for idx in np.split(order, np.flatnonzero(np.diff(codes[order])) + 1):
        k = int(positive[idx].sum())
        if min(k, len(idx) - k) >= 2:
            chunks.append((distance[idx], positive[idx], min(k, len(idx) - k)))
    if not chunks:
        return dict(
            overlap_tiles=0,
            supported_cells=0,
            delta_distance_um=np.nan,
            null_lo=np.nan,
            null_hi=np.nan,
        )
    weights = np.array([c[2] for c in chunks])
    weights = weights / weights.sum()
    observed = sum(
        w * (d[pos].mean() - d[~pos].mean()) for w, (d, pos, _) in zip(weights, chunks)
    )
    draws = np.zeros(199)
    for w, (d, pos, _) in zip(weights, chunks):
        n = len(d)
        k = int(pos.sum())
        total = d.sum()
        for i in range(len(draws)):
            chosen = d[rng.choice(n, k, replace=False)].sum()
            draws[i] += w * (chosen / k - (total - chosen) / (n - k))
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return dict(
        overlap_tiles=len(chunks),
        supported_cells=sum(len(d) for d, _, _ in chunks),
        delta_distance_um=observed,
        null_lo=lo,
        null_hi=hi,
    )


def run():
    start = time.perf_counter()
    genes = pd.read_csv(ROOT / "artifacts/m5/gene_sets/program_genes.csv")
    samples, coverage, rows, inputs = [], [], [], []
    paths = sorted(RAW.glob("*cell_feature_matrix.h5"))
    if len(paths) != 4:
        raise ValueError("Expected the four published matrices")
    for path in paths:
        name = path.name.removesuffix("_cell_feature_matrix.h5")
        coords_path = RAW / f"{name}_cells.parquet.gz"
        with h5py.File(path) as f:
            g = f["matrix"]
            syms = pd.Index(g["features/name"].asstr()[:])
            actual = g["features/feature_type"].asstr()[:] == "Gene Expression"
            measured = set(syms[actual])
            col = syms[actual].get_indexer(PANEL)
            if min(col) < 0 or len(measured) != actual.sum():
                raise ValueError(
                    f"Missing or ambiguous spatial identity markers in {name}"
                )
            # 10x stores genes by cells in CSC; the same arrays encode the transposed CSR.
            X = sparse.csr_matrix(
                (g["data"][:], g["indices"][:], g["indptr"][:]),
                shape=tuple(g["shape"][:][::-1]),
            )
            ids = pd.Index(g["barcodes"].asstr()[:])
            values = pd.DataFrame(
                X[:, np.flatnonzero(actual)[col]].toarray(), columns=PANEL, index=ids
            )
        with gzip.open(coords_path, "rb") as f:
            coords = (
                pq.read_table(f, columns=["cell_id", "x_centroid", "y_centroid"])
                .to_pandas()
                .set_index("cell_id")
            )
        if (
            not ids.is_unique
            or not coords.index.is_unique
            or set(ids) != set(coords.index)
        ):
            raise ValueError("Expression and centroid cell identifiers differ")
        xy = coords.loc[ids].to_numpy()
        if not np.isfinite(xy).all():
            raise ValueError("Nonfinite cell centroid")
        patient = re.search(r"M(\d+)[TN]", name).group(1)
        tissue = "NL" if "AdjNormal" in name else "LM"
        base = dict(sample=name, donor_uid="GSE335552:M" + patient, tissue=tissue)
        for file in (path, coords_path):
            with file.open("rb") as f:
                sha = hashlib.file_digest(f, "sha256").hexdigest()
            inputs.append(
                dict(**base, file=file.name, bytes=file.stat().st_size, sha256=sha)
            )
        for (lineage, program), ss in genes.groupby(["lineage", "program"]):
            frozen = set(ss.gene_symbol)
            covered = sorted(measured & frozen)
            coverage.append(
                dict(
                    **base,
                    lineage=lineage,
                    program=program,
                    n_frozen=len(frozen),
                    n_measured=len(covered),
                    fraction=len(covered) / len(frozen),
                    measured_symbols=" ".join(covered),
                )
            )
        v = values > 0
        t = (v.CD3D | v.CD3E) & (v.CD8A | v.CD8B) & ~v.EPCAM
        my = (
            (v[["CD68", "CD163", "C1QA", "CD14", "LYZ"]].sum(axis=1) >= 2)
            & (v.CD68 | v.CD163 | v.C1QA)
            & ~v.EPCAM
        )
        discordant = t & my
        t = t & ~discordant
        my = my & ~discordant
        samples.append(
            dict(
                **base,
                n_cells=len(v),
                panel_size=len(measured),
                n_myeloid=int(my.sum()),
                n_cd3_cd8=int(t.sum()),
                n_conflicts=int(discordant.sum()),
                n_myeloid_spp1=int((my & v.SPP1).sum()),
                n_myeloid_cxcl9=int((my & v.CXCL9).sum()),
                n_myeloid_cxcl10=int((my & v.CXCL10).sum()),
            )
        )
        if not t.any() or not my.any():
            continue
        dist = cKDTree(xy[t]).query(xy[my], k=1)[0]
        for threshold in (1, 2):
            pos = values.SPP1[my].to_numpy() >= threshold
            for tile in (100, 200, 400):
                stats = spatial_difference(
                    xy[my], dist, pos, tile, np.random.default_rng(20261121)
                )
                rows.append(
                    dict(**base, spp1_min_transcripts=threshold, tile_um=tile, **stats)
                )
    OUT.mkdir(parents=True, exist_ok=True)
    for name, records in [
        ("sample_support", samples),
        ("measured_coverage", coverage),
        ("local_distances", rows),
        ("source_inputs", inputs),
    ]:
        pd.DataFrame(records).to_csv(OUT / f"{name}.csv", index=False)
    print(
        "samples",
        len(samples),
        "patients",
        len(set(s["donor_uid"] for s in samples)),
        "seconds",
        round(time.perf_counter() - start, 2),
    )


if __name__ == "__main__":
    run()
