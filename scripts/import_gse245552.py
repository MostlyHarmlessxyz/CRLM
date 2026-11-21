"""Import the four GSE245552 tissues of one audited donor on CPU."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import resource
import shlex
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from crlm.annotation_m1 import classify

AUDIT = ROOT / "starter/audit_data"
RECORDS = AUDIT / "geo/GSE245552_matrix_audit"
MARKERS = AUDIT / "GSE245552_author_marker_tables.csv"
SEEDS = ROOT / "starter/literature_seeds/literature_program_seeds.csv"
COHORT = "GSE245552"
TISSUES = ("PT", "LM", "NC", "NL")
SEED = 20261121


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def select(manifest, donor):
    with manifest.open(encoding="utf-8-sig", newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r["dataset"] == COHORT and r["donor_uid"] == donor]
    if not rows or len({r["sample_uid"] for r in rows}) != len(rows):
        raise ValueError(f"No samples or duplicate GSMs for {donor}")
    selected = []
    for row in sorted(rows, key=lambda r: (TISSUES.index(r["tissue"]), r["sample_uid"])):
        gsm = row["sample_uid"]
        record = json.loads((RECORDS / f"{gsm}.json").read_text())
        if (record["gsm"] != gsm or record["specimen"] != row["original_specimen_id"]
                or record["tissue"] != row["tissue"] or donor != f"{COHORT}:{record['patient_id']}"
                or row["lesion_uid"] != f"{COHORT}:{record['specimen']}" or row["matrix_url"] != record["matrix_url"]
                or int(row["n_genes_available"]) != record["n_genes"]
                or str(row["patient_mapping_verified"]).lower() != "true" or row["role"] != "development"
                or record["status"] != "verified" or row["input_scale"] != "raw_counts"
                or int(row["input_barcodes"]) != record["n_cells"]):
            raise ValueError(f"{gsm} manifest and source identity or scale mismatch")
        selected.append((row, record))
    return selected


def acquire(record, raw_root):
    directory = raw_root / record["gsm"]
    directory.mkdir(parents=True, exist_ok=True)
    files = {}
    for kind in ("barcode", "feature", "matrix"):
        url = record[f"{kind}_url"]
        path = directory / url.rsplit("/", 1)[-1]
        expected = record[f"{kind}_sha256"]
        if path.exists():
            if sha256(path) != expected:
                raise ValueError(f"Existing input differs from source audit: {path}")
        else:
            with tempfile.NamedTemporaryFile(dir=directory, suffix=".partial") as temp:
                subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location",
                                "--retry", "3", "--retry-all-errors", "--connect-timeout", "15",
                                "--max-time", "240", "--output", temp.name, url], check=True)
                if sha256(Path(temp.name)) != expected:
                    raise ValueError(f"Downloaded file differs from source audit: {url}")
                os.link(temp.name, path)
        files[kind] = path
    return files


def read_library(files, record):
    genes = pd.read_csv(files["feature"], sep="\t", header=None, dtype=str, keep_default_na=False,
                        names=["source_feature_id", "gene_symbol", "feature_type"])
    barcodes = pd.read_csv(files["barcode"], sep="\t", header=None, dtype=str,
                           keep_default_na=False)[0]
    if (genes.eq("").any().any() or genes.source_feature_id.duplicated().any()
            or genes.gene_symbol.duplicated().any() or barcodes.duplicated().any()
            or barcodes.eq("").any() or not genes.feature_type.eq("Gene Expression").all()):
        raise ValueError(f"{record['gsm']} invalid feature or barcode axis")
    with gzip.open(files["matrix"], "rb") as stream, threadpool_limits(limits=4):
        if stream.readline().strip() != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError("Expected a coordinate integer count matrix")
        stream.seek(0)
        coo = mmread(stream, spmatrix=True).tocoo()
    if (coo.shape != (len(genes), len(barcodes)) or coo.shape != (record["n_genes"], record["n_cells"])
            or coo.nnz != record["declared_nnz"] or coo.data.dtype.kind not in "iu"
            or coo.data.min() < 0 or coo.data.max() > np.iinfo(np.int32).max):
        raise ValueError(f"{record['gsm']} invalid matrix shape, entry count or values")
    if coo.row.min() < 0 or coo.row.max() >= len(genes) or coo.col.min() < 0 or coo.col.max() >= len(barcodes):
        raise ValueError("MTX coordinate outside feature or cell axis")
    # Derive reference values from the source COO before transposition.
    rng = np.random.default_rng(SEED)
    pick = rng.choice(coo.nnz, 12, replace=False)
    reference = list(zip(coo.col[pick].tolist(), coo.row[pick].tolist(), coo.data[pick].tolist()))
    for cell, gene in zip(rng.integers(len(barcodes), size=4), rng.integers(len(genes), size=4)):
        match = (coo.col == cell) & (coo.row == gene)
        reference.append((int(cell), int(gene), int(coo.data[match].sum())))
    counts = sparse.csr_matrix(coo.T, dtype=np.int32)
    if counts.nnz != coo.nnz:
        raise ValueError("Duplicate source coordinates")
    counts.eliminate_zeros()
    if int(counts.sum(dtype=np.int64)) != record["total_molecules"]:
        raise ValueError("Visible count total differs from the audited source")
    genes.index = pd.Index(genes.source_feature_id, name="feature_id")
    return counts, genes, barcodes, reference


def qc_cells(counts, symbols):
    totals = np.asarray(counts.sum(axis=1, dtype=np.int64)).ravel()
    detected = np.diff(counts.indptr)
    mt = symbols.str.startswith("MT-").to_numpy()
    pct_mt = np.full(len(totals), np.nan)
    if mt.any():
        mt_counts = np.asarray(counts[:, mt].sum(axis=1, dtype=np.int64)).ravel()
        np.divide(100.0 * mt_counts, totals, out=pct_mt, where=totals > 0)
    frame = pd.DataFrame({"n_counts_visible": totals, "n_genes_visible": detected, "pct_mt": pct_mt,
                          "qc_retained": totals > 0, "qc_exclusion": np.where(totals > 0, "", "empty_library_cell")})
    # Distribution tails are review flags, preserving cells for the coverage denominator.
    for column in ("n_counts_visible", "n_genes_visible"):
        lo, hi = frame[column].quantile([0.01, 0.99])
        frame[f"{column}_tail"] = (frame[column] < lo) | (frame[column] > hi)
    return frame, int(mt.sum())


def write_checked(adata, path, reference):
    with tempfile.TemporaryDirectory(dir=path.parent) as work:
        candidate = Path(work) / path.name
        adata.write_h5ad(candidate, compression="lzf")
        saved = ad.read_h5ad(candidate)
        if (not sparse.isspmatrix_csr(saved.X) or saved.X.dtype != adata.X.dtype
                or saved.shape != adata.shape or saved.X.nnz != adata.X.nnz
                or saved.X.sum(dtype=np.int64) != adata.X.sum(dtype=np.int64)
                or not saved.obs_names.equals(adata.obs_names) or not saved.var_names.equals(adata.var_names)):
            raise ValueError(f"{path.name} count or axis readback mismatch")
        for column in ("sample_uid", "donor_uid", "lesion_uid", "tissue", "barcode", "lineage_m1", "pct_mt"):
            pd.testing.assert_series_equal(saved.obs[column].astype(adata.obs[column].dtype),
                                           adata.obs[column], check_dtype=False, check_categorical=False)
        for cell, gene, value in reference:
            if int(saved.X[cell, gene]) != value:
                raise ValueError(f"{path.name} source entry mismatch at {cell}, {gene}")
        stored_bytes = candidate.stat().st_size
        candidate.replace(path)
    return stored_bytes


def import_sample(row, record, files, processed, source):
    gsm = row["sample_uid"]
    counts, genes, barcodes, reference = read_library(files, record)
    obs, n_mt = qc_cells(counts, genes.gene_symbol)
    labels, detail = classify(counts, genes.gene_symbol.tolist(), MARKERS)
    obs = pd.concat([obs, detail], axis=1)
    obs["lineage_m1"] = labels
    obs["barcode"] = barcodes
    for key in ("dataset", "sample_uid", "donor_uid", "lesion_uid", "tissue", "role", "original_specimen_id"):
        obs[key] = row[key]
    obs["capture_uid"] = gsm
    obs.index = pd.Index(gsm + ":" + barcodes, name="cell_uid")
    qc = {"sample_uid": gsm, "donor_uid": row["donor_uid"], "tissue": row["tissue"],
          "input_cells": len(obs), "retained_cells": int(obs.qc_retained.sum()),
          "excluded_empty": int((~obs.qc_retained).sum()), "measured_mt_genes": n_mt,
          "pct_mt_scope": "visible_axis" if n_mt else "unavailable_no_MT_features",
          "cells_with_positive_mt_counts": int(obs.pct_mt.gt(0).sum()) if n_mt else None,
          "pct_mt_median": obs.pct_mt.median() if n_mt else None}
    for col in ("n_counts_visible", "n_genes_visible"):
        for suffix, quantile in (("p01", .01), ("median", .50), ("p99", .99)):
            qc[f"{col}_{suffix}"] = float(obs[col].quantile(quantile))
        qc[f"{col}_tail_cells"] = int(obs[f"{col}_tail"].sum())
    adata = ad.AnnData(X=counts, obs=obs, var=genes)
    adata.uns["expression"] = {"location": "X", "scale": "raw_counts", "genes_prefiltered": True,
                               "gene_universe_id": row["gene_universe_id"], "source": record["source_record"],
                               "gene_mapping": "exact deposited ID and symbol, no alias substitutions",
                               "prior_antitumor_therapy": "none reported at cohort level (B07)",
                               "collection_window": "paired dates unverified",
                               "capture_protocol": "CD45+ and CD45- magnetic sorting; 2:1 remix reported in B07"}
    adata.uns["provenance"] = dict(source, **{k: record[k] for k in
        ("barcode_sha256", "matrix_sha256", "feature_sha256", "barcode_url", "matrix_url", "feature_url")})
    path = processed / f"{gsm}.h5ad"
    size = write_checked(adata, path, reference)
    record_out = {"sample_uid": gsm, "donor_uid": row["donor_uid"], "tissue": row["tissue"],
                  "n_cells": counts.shape[0], "n_features": counts.shape[1], "nnz": counts.nnz,
                  "visible_counts": int(counts.sum(dtype=np.int64)), "dtype": str(counts.dtype),
                  "raw_bytes": sum(p.stat().st_size for p in files.values()), "h5ad_bytes": size,
                  "h5ad_path": str(path.relative_to(ROOT)), "h5ad_sha256": sha256(path),
                  "source_entries_checked": len(reference), "axes_counts_readback": "passed"}
    levels = ["T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast", "Nonimmune", "Mixed", "Unknown"]
    if not obs.lineage_m1.isin(levels).all():
        raise ValueError(f"{gsm} invalid lineage")
    coverage, review = [], []
    for level in levels:
        group = obs.loc[obs.qc_retained & obs.lineage_m1.eq(level)]
        coverage.append({"donor_uid": row["donor_uid"], "sample_uid": gsm, "tissue": row["tissue"],
                         "lineage": level, "n_cells": len(group), "denominator_retained": qc["retained_cells"]})
        review.append(group.sample(min(3, len(group)), random_state=SEED))
    return record_out, qc, coverage, pd.concat(review), set(genes.gene_symbol)


def run(manifest, donor):
    start = time.perf_counter()
    started = utc_now()
    config = json.loads((ROOT / "configs/project.json").read_text())
    paths = config["paths"]
    raw_root = ROOT / paths["raw"] / COHORT
    processed = ROOT / paths["processed"] / "m1"
    out = ROOT / paths["artifacts"] / "m1"
    processed.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    selected = select(manifest, donor)
    if [r["tissue"] for r, _ in selected] != list(TISSUES):
        raise ValueError(f"M1 quartet requires one audited sample per tissue for {donor}")
    selection = pd.DataFrame([row for row, _ in selected])
    selection.to_csv(out / "sample_selection.csv", index=False)
    source = {"manifest_sha256": sha256(manifest), "marker_table_sha256": sha256(MARKERS),
              "import_code_sha256": sha256(Path(__file__)),
              "annotation_code_sha256": sha256(ROOT / "src/crlm/annotation_m1.py"),
              "recorded_protocol_source": "B07, Cancer Res 84, 1353-1354, doi:10.1158/0008-5472.CAN-23-2123"}
    ad.settings.allow_write_nullable_strings = True
    matrices, metrics, coverage, reviews, inputs, axes = [], [], [], [], [], {}
    for row, record in selected:
        t0 = time.perf_counter()
        files = acquire(record, raw_root)
        for kind, file in files.items():
            inputs.append({"sample_uid": record["gsm"], "kind": kind,
                           "path": str(file.relative_to(ROOT)), "url": record[f"{kind}_url"],
                           "bytes": file.stat().st_size, "verified_sha256": record[f"{kind}_sha256"]})
        matrix, qc, counts, review, axis = import_sample(row, record, files, processed, source)
        matrix["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
        matrices.append(matrix)
        metrics.append(qc)
        coverage.extend(counts)
        reviews.append(review)
        axes[row["tissue"]] = axis
        print(f"{row['sample_uid']} {row['tissue']}: {matrix['n_cells']} cells, readback passed", flush=True)

    common = set.intersection(*axes.values())
    union = set.union(*axes.values())
    expected_common = set((AUDIT / "GSE245552_feature_intersection.txt").read_text().splitlines())
    expected_extra = set((AUDIT / "GSE245552_normal_only_features_353.txt").read_text().splitlines())
    if (axes["PT"] != axes["LM"] or axes["NC"] != axes["NL"] or common != expected_common
            or union - common != expected_extra or axes["PT"] != common or axes["NC"] != union):
        raise ValueError("Feature measurement domains differ from the audited tissue-specific axes")
    features = pd.DataFrame({"gene_symbol": sorted(union)})
    for tissue in TISSUES:
        features[tissue] = features.gene_symbol.isin(axes[tissue])
    features["shared"] = features.gene_symbol.isin(common)
    features.to_csv(out / "gene_coverage.csv", index=False)
    seed_axes = []
    seeds = pd.read_csv(SEEDS)
    for set_id, rows in seeds.groupby("seed_set_id"):
        genes = set(rows.gene_symbol_original)
        for tissue, axis in axes.items():
            absent = sorted(genes - axis)
            seed_axes.append({"seed_set_id": set_id, "tissue": tissue, "n_seed_genes": len(genes),
                              "n_on_axis": len(genes) - len(absent), "absent_exact_symbols": " ".join(absent)})
    pd.DataFrame(seed_axes).to_csv(out / "seed_feature_coverage.csv", index=False)
    pd.DataFrame(inputs).to_csv(out / "source_inputs.csv", index=False)
    pd.DataFrame(matrices).to_csv(out / "matrix_audit.csv", index=False)
    pd.DataFrame(metrics).to_csv(out / "qc_capabilities.csv", index=False)
    pd.DataFrame(coverage).to_csv(out / "lineage_coverage.csv", index=False)
    pd.concat(reviews).to_csv(out / "marker_spotchecks.csv")
    profile = {"started_utc": started, "ended_utc": utc_now(), "cwd": os.getcwd(),
               "command": shlex.join([sys.executable, *sys.argv]), "python": sys.version,
               "packages": {p: metadata.version(p) for p in ("anndata", "numpy", "scipy", "h5py", "pandas", "threadpoolctl")},
               "selection": donor, "seed": SEED, "backend": "CPU", "source": source,
               "seed_table_sha256": sha256(SEEDS),
               "elapsed_seconds": round(time.perf_counter() - start, 3),
               "peak_rss_bytes_linux": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
               "raw_bytes": sum(m["raw_bytes"] for m in matrices),
               "h5ad_bytes": sum(m["h5ad_bytes"] for m in matrices),
               "largest_staged_h5ad_bytes": max(m["h5ad_bytes"] for m in matrices),
               "decompressed_mtx_temp_bytes": 0}
    (out / "resource_profile.json").write_text(json.dumps(profile, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--donor", required=True)
    args = parser.parse_args()
    run(args.manifest, args.donor)
