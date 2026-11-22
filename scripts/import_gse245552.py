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
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
from crlm.annotation_m1 import HIERARCHICAL_SUBTYPES, LEVELS, classify
from crlm.hgnc_alias_map import (
    HGNC_PROVENANCE_PATH,
    HGNC_TSV_PATH,
    SHARED_NUCLEAR_UNIVERSE,
    canonicalize_symbols,
    normalize_gene_axis,
)

AUDIT = ROOT / "starter/audit_data"
RECORDS = AUDIT / "geo/GSE245552_matrix_audit"
MARKERS = AUDIT / "GSE245552_author_marker_tables.csv"
SEEDS = ROOT / "starter/literature_seeds/literature_program_seeds.csv"
COHORT = "GSE245552"
TISSUES = ("PT", "LM", "NC", "NL")
RESOLVED_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast")
CORE_THREE_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma")
MYELOID_FINE_GROUPS = ("Macrophage", "Monocyte")
MT_DEPLETED_FRACTION_THRESHOLD = 0.01
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


def enrich_var_hgnc(genes):
    """Attach deterministic HGNC canonical symbol and shared-universe flags to var."""
    if "canonical_symbol" in genes.columns and "in_shared_nuclear_245_164" in genes.columns:
        return genes
    excluded_353 = set(SHARED_NUCLEAR_UNIVERSE["excluded_353_canonical"])
    tumor_245 = set(SHARED_NUCLEAR_UNIVERSE["245_tumor_canonical_nuclear"])
    shared_164 = set(SHARED_NUCLEAR_UNIVERSE["245_164_shared_nuclear"])
    shared_all = set(SHARED_NUCLEAR_UNIVERSE["245_164_203608_shared_nuclear"])
    mapped = normalize_gene_axis(genes["gene_symbol"].tolist(), excluded_canonical=excluded_353)
    if not mapped["is_unique_measurable"].all() or not mapped["canonical_symbol"].is_unique:
        raise ValueError("GSE245552 feature axis produced non-unique HGNC canonical symbols")
    enriched = genes.copy()
    enriched["canonical_symbol"] = mapped["canonical_symbol"].to_numpy()
    enriched["hgnc_id"] = mapped["hgnc_id"].to_numpy()
    enriched["hgnc_match_type"] = mapped["match_type"].to_numpy()
    enriched["in_excluded_353"] = mapped["in_excluded_353"].to_numpy()
    enriched["in_tumor_canonical_nuclear_245"] = enriched["canonical_symbol"].isin(tumor_245).to_numpy()
    enriched["in_shared_nuclear_245_164"] = enriched["canonical_symbol"].isin(shared_164).to_numpy()
    enriched["in_shared_nuclear_245_164_203608"] = enriched["canonical_symbol"].isin(shared_all).to_numpy()
    return enriched


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
    genes = enrich_var_hgnc(genes)
    return counts, genes, barcodes, reference


def qc_cells(counts, symbols):
    totals = np.asarray(counts.sum(axis=1, dtype=np.int64)).ravel()
    detected = np.diff(counts.indptr)
    mt = pd.Series(symbols).astype(str).str.startswith("MT-").to_numpy()
    pct_mt = np.full(len(totals), np.nan)
    nonempty = int((totals > 0).sum())
    if mt.any():
        mt_counts = np.asarray(counts[:, mt].sum(axis=1, dtype=np.int64)).ravel()
        np.divide(100.0 * mt_counts, totals, out=pct_mt, where=totals > 0)
        pos_mt = int((mt_counts > 0).sum())
        if nonempty > 0 and (pos_mt / nonempty) < MT_DEPLETED_FRACTION_THRESHOLD:
            pct_mt_scope = "visible_axis_depleted_source_anomaly"
            pct_mt_usable = False
        else:
            pct_mt_scope = "visible_axis"
            pct_mt_usable = True
    else:
        pct_mt_scope = "unavailable_no_MT_features"
        pct_mt_usable = False
    frame = pd.DataFrame({
        "n_counts_visible": totals,
        "n_genes_visible": detected,
        "pct_mt": pct_mt,
        "pct_mt_scope": pct_mt_scope,
        "pct_mt_usable_for_qc": pct_mt_usable,
        "qc_retained": totals > 0,
        "qc_exclusion": np.where(totals > 0, "", "empty_library_cell"),
    })
    # Distribution tails are review flags, preserving cells for the M1 coverage denominator.
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
        obs_check_cols = [
            c for c in (
                "sample_uid", "donor_uid", "lesion_uid", "tissue", "barcode",
                "lineage_m1", "subtype_hierarchical", "myeloid_fine",
                "lineage_m2", "subtype_m2", "qc_retained", "pct_mt",
            )
            if c in adata.obs.columns
        ]
        for column in obs_check_cols:
            pd.testing.assert_series_equal(saved.obs[column].astype(adata.obs[column].dtype),
                                           adata.obs[column], check_dtype=False, check_categorical=False)
        for vcol in ("gene_symbol", "canonical_symbol"):
            if vcol in adata.var.columns:
                pd.testing.assert_series_equal(saved.var[vcol].astype(str),
                                               adata.var[vcol].astype(str), check_dtype=False, check_categorical=False)
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

    mt_scope = str(obs["pct_mt_scope"].iloc[0])
    mt_usable = bool(obs["pct_mt_usable_for_qc"].iloc[0])
    pos_mt_cells = int(obs.pct_mt.gt(0).sum()) if n_mt else None
    retained_n = int(obs.qc_retained.sum())
    qc = {
        "sample_uid": gsm,
        "donor_uid": row["donor_uid"],
        "tissue": row["tissue"],
        "input_cells": len(obs),
        "retained_cells": retained_n,
        "excluded_empty": int((~obs.qc_retained).sum()),
        "measured_mt_genes": n_mt,
        "pct_mt_scope": mt_scope,
        "pct_mt_usable_for_qc": mt_usable,
        "cells_with_positive_mt_counts": pos_mt_cells,
        "mt_positive_cell_fraction": (float(pos_mt_cells / retained_n) if (n_mt and retained_n) else None),
        "cells_pct_mt_gt_20": (int(obs.pct_mt.gt(20.0).sum()) if n_mt else None),
        "pct_mt_median": (float(obs.pct_mt.median()) if n_mt else None),
        "pct_mt_max": (float(obs.pct_mt.max()) if (n_mt and obs.pct_mt.notna().any()) else None),
    }
    for col in ("n_counts_visible", "n_genes_visible"):
        for suffix, quantile in (("p01", .01), ("median", .50), ("p99", .99)):
            qc[f"{col}_{suffix}"] = float(obs[col].quantile(quantile))
        qc[f"{col}_tail_cells"] = int(obs[f"{col}_tail"].sum())

    adata = ad.AnnData(X=counts, obs=obs, var=genes)
    adata.uns["expression"] = {
        "location": "X",
        "scale": "raw_counts",
        "genes_prefiltered": True,
        "gene_universe_id": row["gene_universe_id"],
        "source": record["source_record"],
        "gene_mapping": "exact deposited ID/symbol in var.gene_symbol; deterministic HGNC canonical symbol in var.canonical_symbol",
        "pct_mt_scope": mt_scope,
        "pct_mt_usable_for_qc": mt_usable,
        "prior_antitumor_therapy": "none reported at cohort level (B07)",
        "collection_window": "paired dates unverified",
        "capture_protocol": "CD45+ and CD45- magnetic sorting; 2:1 remix reported in B07",
    }
    adata.uns["provenance"] = dict(source, **{k: record[k] for k in
        ("barcode_sha256", "matrix_sha256", "feature_sha256", "barcode_url", "matrix_url", "feature_url")})
    path = processed / f"{gsm}.h5ad"
    size = write_checked(adata, path, reference)
    record_out = {
        "sample_uid": gsm,
        "donor_uid": row["donor_uid"],
        "tissue": row["tissue"],
        "n_cells": counts.shape[0],
        "n_features": counts.shape[1],
        "nnz": counts.nnz,
        "visible_counts": int(counts.sum(dtype=np.int64)),
        "dtype": str(counts.dtype),
        "raw_bytes": sum(p.stat().st_size for p in files.values()),
        "h5ad_bytes": size,
        "h5ad_path": str(path.relative_to(ROOT)),
        "h5ad_sha256": sha256(path),
        "source_entries_checked": len(reference),
        "axes_counts_readback": "passed",
    }
    if not obs.lineage_m1.isin(LEVELS).all():
        raise ValueError(f"{gsm} invalid lineage")

    ret_obs = obs.loc[obs.qc_retained]
    denom_total_imm = int(ret_obs.lineage_m1.isin(RESOLVED_IMMUNE).sum())
    denom_core3_imm = int(ret_obs.lineage_m1.isin(CORE_THREE_IMMUNE).sum())

    def _cov_entry(kind, parent_lin, grp_name, n_grp):
        p_tot = (n_grp / denom_total_imm) if (parent_lin in RESOLVED_IMMUNE and denom_total_imm > 0) else np.nan
        p_c3 = (n_grp / denom_core3_imm) if (parent_lin in CORE_THREE_IMMUNE and denom_core3_imm > 0) else np.nan
        return {
            "donor_uid": row["donor_uid"],
            "sample_uid": gsm,
            "tissue": row["tissue"],
            "level": kind,
            "lineage": parent_lin,
            "group": grp_name,
            "n_cells": int(n_grp),
            "denominator_retained": retained_n,
            "denominator_total_resolved_immune": denom_total_imm,
            "denominator_core_three_immune": denom_core3_imm,
            "prop_total_resolved_immune": p_tot,
            "prop_core_three_immune": p_c3,
        }

    coverage, review = [], []
    for level in LEVELS:
        group = ret_obs.loc[ret_obs.lineage_m1.eq(level)]
        coverage.append(_cov_entry("lineage", level, level, len(group)))
        review.append(group.sample(min(3, len(group)), random_state=SEED))
    for parent_lin in LEVELS:
        for sub_name in HIERARCHICAL_SUBTYPES[parent_lin]:
            n_sub = int((ret_obs.lineage_m1.eq(parent_lin) & ret_obs.subtype_hierarchical.eq(sub_name)).sum())
            coverage.append(_cov_entry("subtype", parent_lin, sub_name, n_sub))
    for mf_name in MYELOID_FINE_GROUPS:
        n_mf = int((ret_obs.lineage_m1.eq("Mono/Mac/DC") & ret_obs.myeloid_fine.eq(mf_name)).sum())
        coverage.append(_cov_entry("myeloid_fine", "Mono/Mac/DC", mf_name, n_mf))

    return record_out, qc, coverage, pd.concat(review), set(genes.gene_symbol), set(genes.canonical_symbol)


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
    source = {
        "manifest_sha256": sha256(manifest),
        "marker_table_sha256": sha256(MARKERS),
        "import_code_sha256": sha256(Path(__file__)),
        "annotation_code_sha256": sha256(ROOT / "src/crlm/annotation_m1.py"),
        "hgnc_code_sha256": sha256(ROOT / "src/crlm/hgnc_alias_map.py"),
        "hgnc_table_sha256": sha256(HGNC_TSV_PATH),
        "hgnc_provenance_sha256": sha256(HGNC_PROVENANCE_PATH),
        "recorded_protocol_source": "B07, Cancer Res 84, 1353-1354, doi:10.1158/0008-5472.CAN-23-2123",
    }
    ad.settings.allow_write_nullable_strings = True
    matrices, metrics, coverage, reviews, inputs, axes, canon_axes = [], [], [], [], [], {}, {}
    for row, record in selected:
        t0 = time.perf_counter()
        files = acquire(record, raw_root)
        for kind, file in files.items():
            inputs.append({
                "sample_uid": record["gsm"],
                "kind": kind,
                "path": str(file.relative_to(ROOT)),
                "url": record[f"{kind}_url"],
                "bytes": file.stat().st_size,
                "verified_sha256": record[f"{kind}_sha256"],
            })
        matrix, qc, counts, review, axis, canon_axis = import_sample(row, record, files, processed, source)
        matrix["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
        matrices.append(matrix)
        metrics.append(qc)
        coverage.extend(counts)
        reviews.append(review)
        axes[row["tissue"]] = axis
        canon_axes[row["tissue"]] = canon_axis
        print(f"{row['sample_uid']} {row['tissue']}: {matrix['n_cells']} cells, readback passed", flush=True)

    common = set.intersection(*axes.values())
    union = set.union(*axes.values())
    expected_common = set((AUDIT / "GSE245552_feature_intersection.txt").read_text().splitlines())
    expected_extra = set((AUDIT / "GSE245552_normal_only_features_353.txt").read_text().splitlines())
    expected_canon_tumor = set(SHARED_NUCLEAR_UNIVERSE["245_tumor_canonical_nuclear"])
    expected_canon_extra = set(SHARED_NUCLEAR_UNIVERSE["excluded_353_canonical"])
    if (axes["PT"] != axes["LM"] or axes["NC"] != axes["NL"] or common != expected_common
            or union - common != expected_extra or axes["PT"] != common or axes["NC"] != union
            or canon_axes["PT"] != expected_canon_tumor
            or (canon_axes["NC"] - canon_axes["PT"]) != expected_canon_extra):
        raise ValueError("Feature measurement domains differ from the audited tissue-specific axes")

    features = pd.DataFrame({"gene_symbol": sorted(union)})
    features_hgnc = enrich_var_hgnc(features)
    for tissue in TISSUES:
        features_hgnc[tissue] = features_hgnc.gene_symbol.isin(axes[tissue])
    features_hgnc["shared"] = features_hgnc.gene_symbol.isin(common)
    col_order = [
        "gene_symbol", "canonical_symbol", "hgnc_id", "hgnc_match_type",
        *TISSUES, "shared",
        "in_excluded_353", "in_tumor_canonical_nuclear_245",
        "in_shared_nuclear_245_164", "in_shared_nuclear_245_164_203608",
    ]
    features_hgnc[col_order].to_csv(out / "gene_coverage.csv", index=False)

    seed_axes = []
    seeds = pd.read_csv(SEEDS)
    for set_id, rows in seeds.groupby("seed_set_id"):
        genes_raw = set(rows.gene_symbol_original)
        genes_canon = set(canonicalize_symbols(rows.gene_symbol_original))
        for tissue, axis in axes.items():
            absent_exact = sorted(genes_raw - axis)
            absent_canon = sorted(genes_canon - canon_axes[tissue])
            seed_axes.append({
                "seed_set_id": set_id,
                "tissue": tissue,
                "n_seed_genes": len(genes_raw),
                "n_on_axis": len(genes_raw) - len(absent_exact),
                "absent_exact_symbols": " ".join(absent_exact),
                "n_seed_genes_canonical": len(genes_canon),
                "n_on_axis_canonical": len(genes_canon) - len(absent_canon),
                "absent_canonical_symbols": " ".join(absent_canon),
            })
    pd.DataFrame(seed_axes).to_csv(out / "seed_feature_coverage.csv", index=False)
    pd.DataFrame(inputs).to_csv(out / "source_inputs.csv", index=False)
    pd.DataFrame(matrices).to_csv(out / "matrix_audit.csv", index=False)
    pd.DataFrame(metrics).to_csv(out / "qc_capabilities.csv", index=False)
    pd.DataFrame(coverage).to_csv(out / "lineage_coverage.csv", index=False)
    pd.concat(reviews).to_csv(out / "marker_spotchecks.csv")
    profile = {
        "started_utc": started,
        "ended_utc": utc_now(),
        "cwd": os.getcwd(),
        "command": shlex.join([sys.executable, *sys.argv]),
        "python": sys.version,
        "packages": {p: metadata.version(p) for p in ("anndata", "numpy", "scipy", "h5py", "pandas", "threadpoolctl")},
        "selection": donor,
        "seed": SEED,
        "backend": "CPU",
        "source": source,
        "seed_table_sha256": sha256(SEEDS),
        "elapsed_seconds": round(time.perf_counter() - start, 3),
        "peak_rss_bytes_linux": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "raw_bytes": sum(m["raw_bytes"] for m in matrices),
        "h5ad_bytes": sum(m["h5ad_bytes"] for m in matrices),
        "largest_staged_h5ad_bytes": max(m["h5ad_bytes"] for m in matrices),
        "decompressed_mtx_temp_bytes": 0,
    }
    (out / "resource_profile.json").write_text(json.dumps(profile, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--donor", required=True)
    args = parser.parse_args()
    run(args.manifest, args.donor)
