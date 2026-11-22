#!/usr/bin/env python3
"""Prospective zero-exposure external validation of frozen `CRLM_V1.1_245_F1` definitions in GSE178318.

Implements:
1. Verified streaming integer UMI import of `GSE178318_{barcodes,genes,matrix}` (Che et al.,
   Cell Discovery 2021; 6 synchronous Stage IV MSS CRC patients, 15 libraries = 140,281 cells).
2. Strict peripheral-blood isolation: documents all 15 GSMs in `coverage/all_15_gsm_barcode_audit.csv`
   and physically excludes the 3 `PBMC` libraries (`GSM5387664`, `GSM5387671`, `GSM5387674`;
   26,950 cells) so ONLY the 12 paired `PT/LM` tumor specimens (6 PT + 6 LM = 113,331 input cells)
   are materialized into `data/processed/m6/GSE178318/{gsm}.h5ad` and paired analyses.
3. Official HGNC + Ensembl-ID (`ENSG`) gene-axis canonicalization (`33,694 -> 33,658` unique
   measurable canonical genes; `13/13` MT- genes; `18,733/18,733 = 100.0%` coverage of
   `245_164_203608_shared_nuclear` and `512/512 = 100.0%` coverage of all M4 PCA HVGs).
4. Zero-modification M2 cell QC (`n_genes >= 300 & n_counts >= 500`, within-lineage 3-MAD lower
   bound, `joint_high` multiplet gate, `pct_mt <= 20.0%`) PLUS sensitivities aligned to Che et al.
   2021 (`n_genes >= 500 & pct_mt <= 15.0%`), `n_genes >= 800` (aligning GSE164522), and `n >= 20/50/100`.
5. Zero-modification `annotation_m1.classify()` two-stage hierarchical annotation (`0.00%`
   unresolved broad-lineage cells; pure `CD8 T cell`, `CD4 T cell`, `NK cell`, `Macrophage`,
   `Monocyte`, `B cell`, `plasma cell` + dual immune composition denominators).
6. Frozen `CRLM_V1.1_245_F1` program scoring (`debiased_ucell_scores` primary + `log1p_cp10k_debiased`,
   `ucell_raw`, `legacy_positive`, and `drop_liver_Hb_DRB345` ambient-liver decontamination),
   `SECONDARY_DT_178_V1_1` Holm-adjusted testing, symmetric `B` vs `Plasma` UPR decomposition,
   and descriptive `3` treatment-naive (`COL07, COL12, COL16`) vs `3` preoperative-chemotherapy
   (`COL15, COL17, COL18`) subgroup stratification.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
import time
from typing import Iterable, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "starter") not in sys.path:
    sys.path.insert(1, str(ROOT / "starter"))

from adjust_registered_tests import adjust_registered  # noqa: E402
from composition_decomposition import decompose_pair  # noqa: E402
from crlm.annotation_m1 import AUXILIARY_MARKERS, PANELS, classify  # noqa: E402
from crlm.hgnc_alias_map import (  # noqa: E402
    align_axis_to_canonical_targets,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
    normalize_gene_axis,
)
from crlm.model_bags import sample_bags, specimen_immune_composition  # noqa: E402
from crlm.program_scores import (  # noqa: E402
    controls_from_dataframe,
    debiased_ucell_scores,
    legacy_positive_rank_scores,
    log1p_cp10k_scores,
    ucell_scores,
)
from paired_effects import calculate_patient_contrasts, summarize_patients  # noqa: E402

DATASET = "GSE178318"
RAW = ROOT / "data/raw" / DATASET
PROCESSED = ROOT / "data/processed/m6" / DATASET
OUT = ROOT / "artifacts/m6" / DATASET
FROZEN = ROOT / "artifacts/m5"
SHARED_UNIVERSE_NAME = "245_164_203608_shared_nuclear"

# Audited 15-library GSM manifest from starter/audit_data/geo/GSE178318_metadata_only.soft
# and Che et al., Cell Discovery 2021 (DOI: 10.1038/s41421-021-00312-y, Table S1 / Methods).
ALL_15_GSM_SPECS = (
    ("GSM5387660", "COL07_CRC", "COL07", "PT", "treatment_naive", "None", "M", 58, "IVA", "MSS", 17000),
    ("GSM5387661", "COL07_LM", "COL07", "LM", "treatment_naive", "None", "M", 58, "IVA", "MSS", 16100),
    ("GSM5387662", "COL12_CRC", "COL12", "PT", "treatment_naive", "None", "M", 65, "IVA", "MSS", 10750),
    ("GSM5387663", "COL12_LM", "COL12", "LM", "treatment_naive", "None", "M", 65, "IVA", "MSS", 3160),
    ("GSM5387664", "COL12_PBMC", "COL12", "PBMC", "treatment_naive", "None", "M", 65, "IVA", "MSS", 12400),
    ("GSM5387665", "COL15_CRC", "COL15", "PT", "preoperative_chemotherapy", "CAPEOX_x3", "F", 35, "IVA", "MSS", 4250),
    ("GSM5387666", "COL15_LM", "COL15", "LM", "preoperative_chemotherapy", "CAPEOX_x3", "F", 35, "IVA", "MSS", 9700),
    ("GSM5387667", "COL16_CRC", "COL16", "PT", "treatment_naive", "None", "M", 71, "IVA", "MSS", 6000),
    ("GSM5387668", "COL16_LM", "COL16", "LM", "treatment_naive", "None", "M", 71, "IVA", "MSS", 8200),
    ("GSM5387669", "COL17_CRC", "COL17", "PT", "preoperative_chemotherapy", "CAPEOX_x4", "F", 52, "IVA", "MSS", 9035),
    ("GSM5387670", "COL17_LM", "COL17", "LM", "preoperative_chemotherapy", "CAPEOX_x4", "F", 52, "IVA", "MSS", 10000),
    ("GSM5387671", "COL17_PBMC", "COL17", "PBMC", "preoperative_chemotherapy", "CAPEOX_x4", "F", 52, "IVA", "MSS", 6700),
    ("GSM5387672", "COL18_CRC", "COL18", "PT", "preoperative_chemotherapy", "FOLFOX_Bevacizumab_x8", "M", 46, "IVA", "MSS", 8700),
    ("GSM5387673", "COL18_LM", "COL18", "LM", "preoperative_chemotherapy", "FOLFOX_Bevacizumab_x8", "M", 46, "IVA", "MSS", 10436),
    ("GSM5387674", "COL18_PBMC", "COL18", "PBMC", "preoperative_chemotherapy", "FOLFOX_Bevacizumab_x8", "M", 46, "IVA", "MSS", 7850),
)
FORBIDDEN_PBMC_GSMS = ("GSM5387664", "GSM5387671", "GSM5387674")

# Exact zero-modification M2 cell QC constants (identical to scripts/process_gse245552.py).
MIN_GENES_PRIMARY = 300
MIN_COUNTS_PRIMARY = 500
MAX_PCT_MT_VISIBLE = 20.0
AUTHOR_B03_MIN_GENES = 500
AUTHOR_B03_MAX_PCT_MT = 15.0
MT_DEPLETED_FRACTION_THRESHOLD = 0.05
LOW_MAD_SIGMA = 3.0
MULTIPLET_MAD_SIGMA = 3.5
MULTIPLET_SOFT_GENES = 6000
MULTIPLET_SOFT_COUNTS = 40000
MULTIPLET_HARD_GENES = 8500
MULTIPLET_HARD_COUNTS = 80000
MIN_GROUP_FOR_WITHIN_MAD = 50

RESOLVED_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast")
CORE_THREE_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma")
PARENT_LINEAGE_MAP = {
    "T/NK": "T/NK",
    "CD4 T cell": "T/NK",
    "CD8 T cell": "T/NK",
    "Double_neg_or_unsplit_T": "T/NK",
    "NK cell": "T/NK",
    "ILC3": "T/NK",
    "Mono/Mac/DC": "Mono/Mac/DC",
    "monocyte/macrophage": "Mono/Mac/DC",
    "Macrophage": "Mono/Mac/DC",
    "Monocyte": "Mono/Mac/DC",
    "conventional DC": "Mono/Mac/DC",
    "plasmacytoid DC": "Mono/Mac/DC",
    "mature DC": "Mono/Mac/DC",
    "B/Plasma": "B/Plasma",
    "B cell": "B/Plasma",
    "plasma cell": "B/Plasma",
    "Mast": "Mast",
    "mast cell": "Mast",
    "Neutrophil": "Neutrophil",
    "neutrophil": "Neutrophil",
}
FOCAL_SUBTYPES_BY_PARENT = {
    "T/NK": ("CD4 T cell", "CD8 T cell", "NK cell"),
    "Mono/Mac/DC": ("Macrophage", "Monocyte", "monocyte/macrophage"),
    "B/Plasma": ("plasma cell", "B cell"),
}
LIVER_RBC = {
    "ALB", "APOA1", "APOA2", "APOC3", "APOH", "TTR",
    "FGA", "FGB", "FGG", "HP", "HPR", "HBA1", "HBA2", "HBB",
}
SCORE_COLUMNS = ["sample_uid", "lineage", "program", "score", "n_cells"]
SEED = 20261121


def sha(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def qc_bounds(
    values: Sequence[float] | np.ndarray,
    low_sigma: float = LOW_MAD_SIGMA,
    high_sigma: float = MULTIPLET_MAD_SIGMA,
    cap_high: float | None = None,
) -> tuple[float, float]:
    logs = np.log1p(np.asarray(values, dtype=np.float64))
    center = float(np.median(logs))
    scale = 1.4826 * float(np.median(np.abs(logs - center)))
    lo = float(np.expm1(center - low_sigma * scale))
    hi = float(np.expm1(center + high_sigma * scale))
    if cap_high is not None:
        hi = min(hi, float(cap_high))
    return lo, hi


def canonicalize_preserving_axis(symbols: Iterable[object], axis_set: set[str]) -> list[str]:
    """Canonicalize gene symbols while preserving any symbol already on `axis_set` verbatim."""
    mapper = get_hgnc_mapper()
    out: list[str] = []
    seen: set[str] = set()
    for raw in symbols:
        if raw is None:
            continue
        s = str(raw).strip()
        if not s:
            continue
        canon = s if s in axis_set else mapper.resolve_symbol(s).canonical_symbol
        if canon and canon not in seen:
            seen.add(canon)
            out.append(canon)
    return out


def registry_requires_debiased(registry: pd.DataFrame) -> bool:
    if "scoring_method" not in registry.columns or registry.empty:
        return True
    methods = set(registry["scoring_method"].astype(str).str.strip())
    return any(m != "positive_midrank_1500" for m in methods)


def verify_freeze_integrity(freeze_dir: Path) -> int:
    """Verify SHA256 integrity of `freeze_dir` files (strictly resolving relative paths)."""
    freeze_resolved = Path(freeze_dir).resolve()
    root_resolved = ROOT.resolve()
    prog_scores_copy = freeze_resolved / "program_scores.py"
    if not prog_scores_copy.is_file():
        raise FileNotFoundError(f"Missing frozen program_scores.py in {freeze_resolved}")
    if sha(prog_scores_copy) != sha(root_resolved / "src/crlm/program_scores.py"):
        raise ValueError("Current program_scores.py differs from the frozen M5 copy")

    checked = 0
    check_all_root = freeze_resolved == (root_resolved / "artifacts/m5").resolve()
    sha_manifest = freeze_resolved / "frozen_inputs_sha256.txt"
    if sha_manifest.is_file():
        for line in sha_manifest.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            digest_hex, rel_path = line.split(None, 1)
            rel_clean = rel_path.strip()
            if rel_clean == "artifacts/m5" or rel_clean.startswith("artifacts/m5/"):
                target = (freeze_resolved / Path(rel_clean).relative_to("artifacts/m5")).resolve()
                should_check = True
            else:
                target = (root_resolved / rel_clean).resolve()
                should_check = check_all_root or target.is_relative_to(freeze_resolved)
            if should_check:
                if not target.is_file():
                    raise FileNotFoundError(f"Frozen M5 file declared in manifest is missing: {rel_clean}")
                if sha(target) != digest_hex:
                    raise ValueError(f"Frozen M5 checksum mismatch for {rel_clean}")
                checked += 1
        if checked == 0:
            raise ValueError(f"Zero freeze_dir files matched in {sha_manifest}")
    return checked


def load_frozen_matched_controls(
    freeze_dir: Path,
    keys: Iterable[tuple[str, str]],
    registry: pd.DataFrame,
    shared_axis: Sequence[str] | None = None,
    *,
    strict_keys: Iterable[tuple[str, str]] | None = None,
) -> dict[tuple[str, str], dict[str, list[str]]]:
    """Load frozen development matched_controls.csv from `freeze_dir` without external refitting."""
    freeze_resolved = Path(freeze_dir).resolve()
    root_resolved = ROOT.resolve()

    cand_paths: list[Path] = [
        freeze_resolved / "gene_sets/matched_controls.csv",
        freeze_resolved / "matched_controls.csv",
    ]
    if "matched_controls_file" in registry.columns:
        for rel in registry["matched_controls_file"].dropna().unique():
            rel_s = str(rel).strip()
            if not rel_s:
                continue
            if rel_s.startswith("artifacts/m5/"):
                cand_paths.append(freeze_resolved / Path(rel_s).relative_to("artifacts/m5"))
            else:
                cand_paths.append((root_resolved / rel_s).resolve())

    mc_path = next((p for p in cand_paths if p.is_file()), None)
    requires_debiased = registry_requires_debiased(registry)
    if mc_path is None:
        if requires_debiased:
            raise FileNotFoundError(
                f"V1.1 debiased external scoring requires frozen matched_controls.csv under {freeze_resolved}; "
                "refitting control bins on external cohorts is prohibited"
            )
        return {}

    axis_set = (
        set(str(g) for g in shared_axis)
        if shared_axis is not None
        else set(load_shared_nuclear_universe(SHARED_UNIVERSE_NAME))
    )
    mc_df = pd.read_csv(mc_path, dtype=str, keep_default_na=False)
    mc_df["replicate"] = mc_df["replicate"].astype(int)
    mapper = get_hgnc_mapper()
    uniq_genes = set(mc_df["target_gene"]).union(mc_df["control_gene"])
    canon_lookup = {
        g: (g if g in axis_set else mapper.resolve_symbol(g).canonical_symbol)
        for g in uniq_genes
    }
    mc_df["target_gene"] = mc_df["target_gene"].map(canon_lookup)
    mc_df["control_gene"] = mc_df["control_gene"].map(canon_lookup)

    key_list = list(keys)
    resolved: dict[tuple[str, str], dict[str, list[str]]] = {}
    has_lin = "lineage" in mc_df.columns
    for lin, prog in key_list:
        if has_lin:
            sub = mc_df.loc[mc_df["lineage"].eq(lin) & mc_df["program"].eq(prog)]
            if sub.empty and lin in PARENT_LINEAGE_MAP:
                sub = mc_df.loc[mc_df["lineage"].eq(PARENT_LINEAGE_MAP[lin]) & mc_df["program"].eq(prog)]
        else:
            sub = mc_df.loc[mc_df["program"].eq(prog)]
        if not sub.empty:
            parsed = controls_from_dataframe(sub[["program", "target_gene", "replicate", "control_gene"]])
            if prog in parsed and parsed[prog]:
                resolved[(lin, prog)] = parsed[prog]

    if requires_debiased:
        must_have = list(strict_keys) if strict_keys is not None else key_list
        missing_ctrls = [f"{lin}::{prog}" for lin, prog in must_have if (lin, prog) not in resolved]
        if missing_ctrls:
            raise ValueError(
                f"Missing frozen matched_controls in {mc_path} for required programs {missing_ctrls}; "
                "silent fallback to non-debiased scores is prohibited"
            )
    return resolved


def parse_and_audit_barcodes(barcodes: Sequence[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Map every deposited barcode to its verified GSM, patient, tissue, and PBMC exclusion flag."""
    if len(barcodes) == 0 or len(set(barcodes)) != len(barcodes):
        raise ValueError("Deposited GSE178318 barcodes must be non-empty and globally unique")

    spec_by_tag = {
        tag: {
            "gsm": gsm,
            "barcode_tag": tag,
            "patient_id": pat,
            "donor_uid": f"{DATASET}:{pat}",
            "sample_uid": gsm,
            "lesion_uid": f"{DATASET}:{tag}",
            "tissue": tissue,
            "treatment_status": tx,
            "preop_regimen": regimen,
            "sex": sex,
            "age": age,
            "stage": stage,
            "msi_status": msi,
            "expected_cells": expected_n,
            "is_tumor_pt_lm": tissue in ("PT", "LM"),
            "role": "external" if tissue in ("PT", "LM") else "excluded_PBMC",
            "exclusion_reason": "" if tissue in ("PT", "LM") else "peripheral_blood_PBMC_excluded_from_tumor_PT_LM_pairing",
        }
        for (gsm, tag, pat, tissue, tx, regimen, sex, age, stage, msi, expected_n) in ALL_15_GSM_SPECS
    }

    tags = []
    short_bcs = []
    for bc in barcodes:
        if "_" not in bc:
            raise ValueError(f"Unexpected GSE178318 barcode without sample suffix: {bc}")
        short_bc, tag = bc.split("_", 1)
        if tag not in spec_by_tag:
            raise ValueError(f"Unrecognized GSE178318 sample tag {tag!r} in barcode {bc!r}")
        short_bcs.append(short_bc)
        tags.append(tag)

    cell_df = pd.DataFrame(
        {
            "raw_column_index": np.arange(len(barcodes), dtype=np.int32),
            "deposited_barcode": list(barcodes),
            "barcode_16bp": short_bcs,
            "barcode_tag": tags,
        }
    )
    tag_counts = cell_df["barcode_tag"].value_counts().to_dict()

    audit_rows = []
    for (_, tag, *_) in ALL_15_GSM_SPECS:
        info = dict(spec_by_tag[tag])
        observed_n = int(tag_counts.get(tag, 0))
        if observed_n != info["expected_cells"]:
            raise ValueError(
                f"GSE178318 barcode count mismatch for {tag} ({info['gsm']}): "
                f"observed {observed_n} vs expected {info['expected_cells']}"
            )
        info["observed_cells"] = observed_n
        info["eligible_for_paired_analysis"] = bool(info["is_tumor_pt_lm"])
        info["h5ad_materialized"] = bool(info["is_tumor_pt_lm"])
        audit_rows.append(info)

    gsm_audit_df = pd.DataFrame(audit_rows)
    cell_df = cell_df.merge(
        gsm_audit_df[
            [
                "barcode_tag",
                "gsm",
                "sample_uid",
                "donor_uid",
                "lesion_uid",
                "patient_id",
                "tissue",
                "treatment_status",
                "preop_regimen",
                "role",
                "is_tumor_pt_lm",
            ]
        ],
        on="barcode_tag",
        how="left",
        validate="many_to_one",
    ).sort_values("raw_column_index").reset_index(drop=True)
    return cell_df, gsm_audit_df


def stream_tumor_count_matrix(
    mtx_path: Path,
    n_genes: int,
    n_all_cells: int,
    tumor_cell_map: np.ndarray,
) -> tuple[sparse.csr_matrix, int, int]:
    """Stream MatrixMarket integer triplets, verifying full file integrity while retaining only tumor cells."""
    n_tumor = int((tumor_cell_map >= 0).sum())
    blocks = []
    checked = 0
    total_all_umi = 0
    with gzip.open(mtx_path, "rb") as handle:
        banner = handle.readline().strip()
        if banner != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"Expected coordinate integer general MatrixMarket, got {banner!r}")
        header = handle.readline()
        while header.startswith(b"%"):
            header = handle.readline()
        nf, nc, nnz = map(int, header.split())
        if (nf, nc) != (n_genes, n_all_cells):
            raise ValueError(f"MatrixMarket shape {(nf, nc)} disagrees with {(n_genes, n_all_cells)}")

        for chunk in pd.read_csv(
            handle,
            sep=r"\s+",
            header=None,
            names=["g", "c", "count"],
            dtype=np.int64,
            chunksize=2_000_000,
        ):
            g = chunk["g"].to_numpy() - 1
            c = chunk["c"].to_numpy() - 1
            val = chunk["count"].to_numpy()
            if (
                g.min() < 0
                or g.max() >= nf
                or c.min() < 0
                or c.max() >= nc
                or val.min() < 0
                or val.max() > np.iinfo(np.int32).max
            ):
                raise ValueError("Invalid GSE178318 MatrixMarket coordinate or negative/overflow integer count")
            checked += len(val)
            total_all_umi += int(val.sum())
            keep = (tumor_cell_map[c] >= 0) & (val > 0)
            if keep.any():
                blocks.append(
                    (
                        tumor_cell_map[c[keep]],
                        g[keep].astype(np.int32),
                        val[keep].astype(np.int32),
                    )
                )
    if checked != nnz:
        raise ValueError(f"Read {checked} MatrixMarket entries, expected {nnz}")

    rows_arr, cols_arr, vals_arr = (np.concatenate(part) for part in zip(*blocks))
    tumor_csr = sparse.coo_matrix(
        (vals_arr, (rows_arr, cols_arr)),
        shape=(n_tumor, n_genes),
        dtype=np.int32,
    ).tocsr()
    if tumor_csr.nnz != len(vals_arr):
        raise ValueError("Duplicate (cell, gene) coordinates detected in GSE178318 MatrixMarket")
    tumor_csr.eliminate_zeros()
    return tumor_csr, nnz, total_all_umi


def annotate_and_qc_sample(
    X_raw_sample: sparse.csr_matrix,
    sample_meta: Mapping[str, object],
    sample_cells: pd.DataFrame,
    raw_mt_mask: np.ndarray,
    measurable_mask: np.ndarray,
    var_measurable: pd.DataFrame,
) -> tuple[ad.AnnData, dict[str, object]]:
    """Run zero-modification M2 QC and `annotation_m1.classify()` on a single tumor library."""
    gsm = str(sample_meta["sample_uid"])
    n_raw_vars = X_raw_sample.shape[1]
    totals = np.asarray(X_raw_sample.sum(axis=1, dtype=np.int64)).ravel()
    detected = np.diff(X_raw_sample.indptr)

    pct_mt = np.full(len(totals), np.nan, dtype=np.float64)
    nonempty_n = int((totals > 0).sum())
    n_mt_genes = int(raw_mt_mask.sum())
    if n_mt_genes > 0:
        mt_counts = np.asarray(X_raw_sample[:, raw_mt_mask].sum(axis=1, dtype=np.int64)).ravel()
        np.divide(100.0 * mt_counts, totals, out=pct_mt, where=totals > 0)
        pos_mt = int((mt_counts > 0).sum())
        if nonempty_n > 0 and (pos_mt / nonempty_n) < MT_DEPLETED_FRACTION_THRESHOLD:
            pct_mt_scope = "visible_axis_depleted_source_anomaly"
            pct_mt_usable = False
        else:
            pct_mt_scope = "visible_axis"
            pct_mt_usable = True
    else:
        pct_mt_scope = "unavailable_no_MT_features"
        pct_mt_usable = False

    X_meas = X_raw_sample[:, measurable_mask].tocsr()
    X_meas.eliminate_zeros()
    dropped_unmeasurable_umi = int(totals.sum() - int(X_meas.sum(dtype=np.int64)))

    canon_syms = var_measurable["canonical_symbol"].astype(str).tolist()
    labels, detail = classify(X_meas, canon_syms)

    obs = pd.DataFrame(
        {
            "dataset": DATASET,
            "sample_uid": gsm,
            "donor_uid": str(sample_meta["donor_uid"]),
            "lesion_uid": str(sample_meta["lesion_uid"]),
            "patient_id": str(sample_meta["patient_id"]),
            "tissue": str(sample_meta["tissue"]),
            "role": "external",
            "treatment_status": str(sample_meta["treatment_status"]),
            "preop_regimen": str(sample_meta["preop_regimen"]),
            "barcode": sample_cells["deposited_barcode"].to_numpy(dtype=str),
            "barcode_16bp": sample_cells["barcode_16bp"].to_numpy(dtype=str),
            "n_counts_visible": totals,
            "n_genes_visible": detected,
            "pct_mt": pct_mt,
            "pct_mt_scope": pct_mt_scope,
            "pct_mt_usable_for_qc": pct_mt_usable,
        }
    )
    obs = pd.concat([obs, detail.reset_index(drop=True)], axis=1)
    obs["lineage_m1"] = np.asarray(labels, dtype=str)
    obs["lineage_m2"] = obs["lineage_m1"].astype(str)
    obs["subtype_m2"] = obs["subtype_hierarchical"].astype(str)
    obs["myeloid_fine"] = obs["myeloid_fine"].fillna("").astype(str)
    if obs["subtype_m2"].str.contains("unresolved", case=False, regex=False).any():
        raise ValueError(f"{gsm} contains unresolved subtype labels after hierarchical disambiguation")

    low_counts, high_counts = qc_bounds(obs["n_counts_visible"])
    low_genes, high_genes = qc_bounds(obs["n_genes_visible"], cap_high=n_raw_vars)
    obs["qc_joint_low"] = (obs["n_counts_visible"] < low_counts) & (obs["n_genes_visible"] < low_genes)

    within_low = pd.Series(False, index=obs.index)
    within_high = pd.Series(False, index=obs.index)
    for _, grp in obs.groupby("lineage_m2", observed=True):
        if len(grp) >= MIN_GROUP_FOR_WITHIN_MAD:
            lo_c, hi_c = qc_bounds(grp["n_counts_visible"])
            lo_g, hi_g = qc_bounds(grp["n_genes_visible"], cap_high=n_raw_vars)
        else:
            lo_c, hi_c = low_counts, high_counts
            lo_g, hi_g = low_genes, high_genes
        within_low.loc[grp.index] = (grp["n_counts_visible"] < lo_c) & (grp["n_genes_visible"] < lo_g)
        within_high.loc[grp.index] = (
            (grp["n_counts_visible"] > hi_c)
            & (grp["n_genes_visible"] > hi_g)
            & ((grp["n_genes_visible"] > MULTIPLET_SOFT_GENES) | (grp["n_counts_visible"] > MULTIPLET_SOFT_COUNTS))
        )

    hard_high = (obs["n_genes_visible"] > MULTIPLET_HARD_GENES) | (obs["n_counts_visible"] > MULTIPLET_HARD_COUNTS)
    obs["qc_within_lineage_low"] = within_low.to_numpy(dtype=bool)
    obs["qc_joint_high"] = (within_high | hard_high).to_numpy(dtype=bool)

    empty = obs["n_counts_visible"].eq(0).to_numpy(dtype=bool)
    below_floor = (
        (obs["n_genes_visible"] < MIN_GENES_PRIMARY) | (obs["n_counts_visible"] < MIN_COUNTS_PRIMARY)
    ).to_numpy(dtype=bool)
    obs["qc_below_min_floor"] = below_floor

    high_mt = (obs["pct_mt"].to_numpy() > MAX_PCT_MT_VISIBLE) if pct_mt_usable else np.zeros(len(obs), dtype=bool)
    obs["qc_high_mt"] = high_mt

    excluded = empty | below_floor | obs["qc_within_lineage_low"].to_numpy() | obs["qc_joint_high"].to_numpy() | high_mt
    obs["qc_retained"] = ~excluded
    obs["qc_exclusion"] = np.select(
        [empty, below_floor, obs["qc_within_lineage_low"].to_numpy(), obs["qc_joint_high"].to_numpy(), high_mt],
        ["empty", "below_min_300g_500c", "within_lineage_low", "joint_high_multiplet", "high_pct_mt_gt20"],
        default="",
    )

    w_low_arr = obs["qc_within_lineage_low"].to_numpy()
    j_high_arr = obs["qc_joint_high"].to_numpy()
    g_arr = obs["n_genes_visible"].to_numpy()
    c_arr = obs["n_counts_visible"].to_numpy()
    mt_arr = obs["pct_mt"].to_numpy()

    obs["qc_pass_g200"] = ~(empty | (g_arr < 200) | (c_arr < 400) | w_low_arr | j_high_arr | high_mt)
    obs["qc_pass_g300_primary"] = obs["qc_retained"].to_numpy(dtype=bool)
    obs["qc_pass_g500"] = ~(empty | (g_arr < 500) | (c_arr < 500) | w_low_arr | j_high_arr | high_mt)
    obs["qc_pass_g800"] = ~(empty | (g_arr < 800) | (c_arr < 800) | w_low_arr | j_high_arr | high_mt)
    obs["qc_author_b03_500g_15mt"] = ~(
        empty
        | (g_arr < AUTHOR_B03_MIN_GENES)
        | (c_arr < MIN_COUNTS_PRIMARY)
        | w_low_arr
        | j_high_arr
        | (mt_arr > AUTHOR_B03_MAX_PCT_MT if pct_mt_usable else False)
    )
    obs["qc_pass_no_mt_filter"] = ~(empty | below_floor | w_low_arr | j_high_arr)

    obs.index = pd.Index(
        gsm + ":" + sample_cells["deposited_barcode"].astype(str).to_numpy(),
        name="cell_uid",
    )
    if not obs.index.is_unique:
        raise ValueError(f"Duplicate cell_uid in {gsm}")

    adata = ad.AnnData(X=X_meas, obs=obs, var=var_measurable.copy())
    adata.uns["qc_limits"] = {
        "low_counts": low_counts,
        "high_counts": high_counts,
        "low_genes": low_genes,
        "high_genes": high_genes,
        "min_genes_primary": MIN_GENES_PRIMARY,
        "min_counts_primary": MIN_COUNTS_PRIMARY,
        "max_pct_mt_visible": MAX_PCT_MT_VISIBLE,
        "author_b03_min_genes": AUTHOR_B03_MIN_GENES,
        "author_b03_max_pct_mt": AUTHOR_B03_MAX_PCT_MT,
        "pct_mt_scope": pct_mt_scope,
        "pct_mt_usable_for_qc": pct_mt_usable,
    }
    adata.uns["expression"] = {
        "location": "X",
        "scale": "raw_counts",
        "raw_deposited_genes": int(n_raw_vars),
        "unique_measurable_canonical_genes": int(adata.n_vars),
        "shared_nuclear_universe_covered": SHARED_UNIVERSE_NAME,
        "capture_protocol": "Che et al. 2021 (B03) fresh unsorted viable single cells (Calcein AM+ / TO-PRO-3-), 10x 3' v2 GRCh38",
        "treatment_status": str(sample_meta["treatment_status"]),
        "preop_regimen": str(sample_meta["preop_regimen"]),
    }

    keep = obs["qc_retained"].to_numpy(dtype=bool)
    broad = obs["lineage_m2"].to_numpy(dtype=str)
    qc_summary = {
        "sample_uid": gsm,
        "donor_uid": str(sample_meta["donor_uid"]),
        "lesion_uid": str(sample_meta["lesion_uid"]),
        "patient_id": str(sample_meta["patient_id"]),
        "tissue": str(sample_meta["tissue"]),
        "treatment_status": str(sample_meta["treatment_status"]),
        "preop_regimen": str(sample_meta["preop_regimen"]),
        "input_cells": int(len(obs)),
        "retained_cells_m2_primary": int(keep.sum()),
        "retained_cells_author_b03_500g_15mt": int(obs["qc_author_b03_500g_15mt"].sum()),
        "retained_cells_g500": int(obs["qc_pass_g500"].sum()),
        "retained_cells_g800_align164": int(obs["qc_pass_g800"].sum()),
        "excluded_total": int((~keep).sum()),
        "excluded_empty": int((obs["qc_exclusion"] == "empty").sum()),
        "excluded_below_min_300g_500c": int((obs["qc_exclusion"] == "below_min_300g_500c").sum()),
        "excluded_within_lineage_low": int((obs["qc_exclusion"] == "within_lineage_low").sum()),
        "excluded_joint_high_multiplet": int((obs["qc_exclusion"] == "joint_high_multiplet").sum()),
        "excluded_high_pct_mt_gt20": int((obs["qc_exclusion"] == "high_pct_mt_gt20").sum()),
        "n_total_resolved_immune": int((keep & np.isin(broad, RESOLVED_IMMUNE)).sum()),
        "n_core_three_lineage_immune": int((keep & np.isin(broad, CORE_THREE_IMMUNE)).sum()),
        "nonimmune_retained": int((keep & (broad == "Nonimmune")).sum()),
        "mixed_retained": int((keep & (broad == "Mixed")).sum()),
        "unknown_retained": int((keep & (broad == "Unknown")).sum()),
        "raw_total_umi": int(totals.sum()),
        "measurable_total_umi": int(X_meas.sum(dtype=np.int64)),
        "unmeasurable_duplicate_symbol_umi": dropped_unmeasurable_umi,
        "measured_mt_genes": n_mt_genes,
        "pct_mt_median_all": float(np.nanmedian(pct_mt)),
        "pct_mt_median_retained": float(obs.loc[keep, "pct_mt"].median()) if keep.any() else np.nan,
        "n_genes_median_retained": float(obs.loc[keep, "n_genes_visible"].median()) if keep.any() else np.nan,
        "n_counts_median_retained": float(obs.loc[keep, "n_counts_visible"].median()) if keep.any() else np.nan,
    }
    return adata, qc_summary


def write_checked_h5ad(adata: ad.AnnData, path: Path) -> int:
    """Atomically write sparse `.h5ad` and verify CSR integer readback & unique canonical gene axis."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent) as work:
        cand = Path(work) / path.name
        adata.write_h5ad(cand, compression="lzf")
        saved = ad.read_h5ad(cand)
        if (
            not sparse.isspmatrix_csr(saved.X)
            or saved.X.dtype != adata.X.dtype
            or saved.shape != adata.shape
            or saved.X.nnz != adata.X.nnz
            or int(saved.X.sum(dtype=np.int64)) != int(adata.X.sum(dtype=np.int64))
            or not saved.obs_names.equals(adata.obs_names)
            or not saved.var_names.equals(adata.var_names)
            or not saved.var["canonical_symbol"].is_unique
        ):
            raise ValueError(f"Readback integrity check failed for {path.name}")
        size_bytes = cand.stat().st_size
        cand.replace(path)
    return size_bytes


def leave_one_out(patients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (lineage, program), group in patients.groupby(["lineage", "program"]):
        available = group.dropna(subset=["D_T"])
        for patient in available["donor_uid"]:
            others = available.loc[available["donor_uid"] != patient, "D_T"]
            rows.append(
                (
                    lineage,
                    program,
                    patient,
                    len(others),
                    float(others.mean()) if len(others) else np.nan,
                )
            )
    return pd.DataFrame(
        rows, columns=["lineage", "program", "omitted_patient", "n_remaining", "mean_D_T"]
    )


def decompose_subtypes(manifest: pd.DataFrame, subtype_sample_df: pd.DataFrame) -> pd.DataFrame:
    """Lesion-equal-weighted symmetric composition/state decomposition via starter.decompose_pair."""
    usable = subtype_sample_df.loc[subtype_sample_df["n_cells"] >= 50].merge(
        manifest[["sample_uid", "donor_uid", "lesion_uid", "tissue"]],
        on="sample_uid",
        validate="many_to_one",
    )
    usable = usable.loc[usable["tissue"].isin(["PT", "LM"])].copy()
    if usable.empty:
        return pd.DataFrame(
            columns=[
                "donor_uid", "lineage", "program", "min_subtype_cells", "status",
                "observed_mean_PT", "observed_mean_LM", "observed_total_difference",
                "full_composition", "full_within_type_state",
                "common_retained_mass_PT", "common_retained_mass_LM",
                "common_composition", "common_within_type_state", "unidentified_subtypes",
            ]
        )
    usable["proportion"] = usable["n_subtype"] / usable["n_cells"]
    usable["score_mass"] = usable["proportion"] * usable["mean_score"]
    rows = []
    for (donor_uid, lineage, program), group in usable.groupby(["donor_uid", "lineage", "program"]):
        frame = []
        observed: dict[str, float] = {}
        names = sorted(group["subtype"].unique())
        for tissue, tissue_rows in group.groupby("tissue"):
            n_lesions = tissue_rows["lesion_uid"].nunique()
            agg = (
                tissue_rows.groupby("subtype")[["n_subtype", "proportion", "score_mass"]]
                .sum()
                .reindex(names, fill_value=0)
            )
            p = agg["proportion"] / n_lesions
            mu = (agg["score_mass"] / agg["proportion"]).where(p > 0)
            observed[tissue] = float(agg["score_mass"].sum() / n_lesions)
            frame.extend(
                {
                    "tissue": tissue,
                    "subtype": s,
                    "proportion": float(p.at[s]),
                    "mean_score": float(mu.at[s]) if pd.notna(mu.at[s]) else np.nan,
                    "n_subtype": int(agg["n_subtype"].at[s]),
                }
                for s in names
            )
        values = pd.DataFrame(frame)
        unmasked_pt = observed.get("PT")
        unmasked_lm = observed.get("LM")
        unmasked_diff = (
            (unmasked_lm - unmasked_pt)
            if (unmasked_pt is not None and unmasked_lm is not None)
            else None
        )
        for support_min in (1, 20):
            checked = values.copy()
            checked.loc[checked["n_subtype"] < support_min, "mean_score"] = np.nan
            pair = decompose_pair(checked)
            common = pair.get("common_support") or {}
            full = pair.get("full_decomposition") or {}
            cond = common.get("renormalized_components") or {}
            rows.append(
                {
                    "donor_uid": donor_uid,
                    "lineage": lineage,
                    "program": program,
                    "min_subtype_cells": support_min,
                    "status": pair["status"],
                    "observed_mean_PT": unmasked_pt if pair["status"] != "missing_tissue" else None,
                    "observed_mean_LM": unmasked_lm if pair["status"] != "missing_tissue" else None,
                    "observed_total_difference": unmasked_diff if pair["status"] != "missing_tissue" else None,
                    "full_composition": full.get("composition"),
                    "full_within_type_state": full.get("within_type_state"),
                    "common_retained_mass_PT": common.get("retained_mass_PT"),
                    "common_retained_mass_LM": common.get("retained_mass_LM"),
                    "common_composition": cond.get("composition"),
                    "common_within_type_state": cond.get("within_type_state"),
                    "unidentified_subtypes": "|".join(pair.get("unidentified_subtypes") or []),
                }
            )
    return pd.DataFrame(rows)


def evaluate_tumor_cohort(
    X_shared: sparse.csr_matrix,
    shared_axis: Sequence[str],
    obs_all: pd.DataFrame,
    samples: pd.DataFrame,
    eval_registry: pd.DataFrame,
    sets: pd.DataFrame,
    controls_map: Mapping[tuple[str, str], Mapping[str, Sequence[str]]],
    *,
    requires_debiased: bool = True,
) -> dict[str, pd.DataFrame]:
    """Evaluate frozen V1.1 programs, sensitivities, focal subtypes, decompositions, and treatment subgroups."""
    shared_list = [str(g) for g in shared_axis]
    shared_set = set(shared_list)
    shared_idx = pd.Index(shared_list)

    broad = obs_all["lineage_m2"].to_numpy(dtype=str)
    sub = obs_all["subtype_m2"].to_numpy(dtype=str)
    mf = obs_all["myeloid_fine"].to_numpy(dtype=str)
    sample_arr = obs_all["sample_uid"].to_numpy(dtype=str)

    m_primary = obs_all["qc_retained"].to_numpy(dtype=bool)
    m_b03 = obs_all["qc_author_b03_500g_15mt"].to_numpy(dtype=bool)
    m_g500 = obs_all["qc_pass_g500"].to_numpy(dtype=bool)
    m_g800 = obs_all["qc_pass_g800"].to_numpy(dtype=bool)
    m_g200 = obs_all["qc_pass_g200"].to_numpy(dtype=bool)
    m_nomt = obs_all["qc_pass_no_mt_filter"].to_numpy(dtype=bool)
    any_qc = m_primary | m_b03 | m_g500 | m_g800 | m_g200 | m_nomt

    # Dual-denominator sample immune composition & lineage support tables
    comp_rows = []
    support_rows = []
    for s in samples.itertuples(index=False):
        sm = (sample_arr == s.sample_uid) & m_primary
        n_tot_imm = int((sm & np.isin(broad, RESOLVED_IMMUNE)).sum())
        n_c3_imm = int((sm & np.isin(broad, CORE_THREE_IMMUNE)).sum())
        n_ret = int(sm.sum())
        for lin_name in ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast", "Nonimmune", "Mixed", "Unknown"):
            n_lin = int((sm & (broad == lin_name)).sum())
            comp_rows.append(
                {
                    "sample_uid": s.sample_uid,
                    "donor_uid": s.donor_uid,
                    "tissue": s.tissue,
                    "treatment_status": s.treatment_status,
                    "lineage": lin_name,
                    "n_cells_retained": n_lin,
                    "n_sample_retained": n_ret,
                    "n_total_resolved_immune": n_tot_imm,
                    "n_core_three_immune": n_c3_imm,
                    "prop_sample_retained": (n_lin / n_ret if n_ret > 0 else np.nan),
                    "prop_total_resolved_immune": (
                        n_lin / n_tot_imm if n_tot_imm > 0 and lin_name in RESOLVED_IMMUNE else np.nan
                    ),
                    "prop_core_three_immune": (
                        n_lin / n_c3_imm if n_c3_imm > 0 and lin_name in CORE_THREE_IMMUNE else np.nan
                    ),
                }
            )
            if lin_name in CORE_THREE_IMMUNE:
                support_rows.append(
                    {
                        "sample_uid": s.sample_uid,
                        "donor_uid": s.donor_uid,
                        "tissue": s.tissue,
                        "treatment_status": s.treatment_status,
                        "lineage": lin_name,
                        "n_cells_m2_primary": n_lin,
                        "n_cells_author_b03_500g_15mt": int(((sample_arr == s.sample_uid) & m_b03 & (broad == lin_name)).sum()),
                        "n_cells_g800_align164": int(((sample_arr == s.sample_uid) & m_g800 & (broad == lin_name)).sum()),
                        "eligible_n20": n_lin >= 20,
                        "eligible_n50": n_lin >= 50,
                        "eligible_n100": n_lin >= 100,
                    }
                )

    measured_rows = []
    cell_score_frames = []

    for r in eval_registry.itertuples(index=False):
        lin = str(r.lineage)
        prog = str(r.program)
        raw_req = sets.loc[sets["program"].eq(prog) & sets["lineage"].eq(lin), "gene_symbol"].tolist()
        canon_req = canonicalize_preserving_axis(raw_req, shared_set)
        keep = [g for g in canon_req if g in shared_set]
        ok = (
            len(keep) >= int(r.min_measured_positive_genes)
            and len(keep) >= len(canon_req) * float(r.min_measured_fraction_positive)
        )

        if lin in ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast", "Neutrophil"):
            in_lin = broad == lin
        else:
            parent = PARENT_LINEAGE_MAP.get(lin, lin)
            in_lin = (broad == parent) & ((sub == lin) | (mf == lin))

        measured_rows.append(
            {
                "lineage": lin,
                "program": prog,
                "n_axis_genes": len(shared_list),
                "n_frozen_raw": len(raw_req),
                "n_frozen_canonical": len(canon_req),
                "n_measured": len(keep),
                "measured_fraction": (len(keep) / len(canon_req) if canon_req else 0.0),
                "has_frozen_matched_controls": (lin, prog) in controls_map,
                "n_cells_m2_primary": int((m_primary & in_lin).sum()),
                "n_cells_author_b03_500g_15mt": int((m_b03 & in_lin).sum()),
                "missing_symbols": "|".join(sorted(set(canon_req) - set(keep))),
                "eligible": ok,
            }
        )

        if ok and requires_debiased and (lin, prog) not in controls_map:
            raise ValueError(
                f"Missing frozen matched_controls for {lin}::{prog}; silent fallback to non-debiased scores is prohibited"
            )

        union_rows = np.flatnonzero(any_qc & in_lin)
        if not ok or len(union_rows) == 0:
            continue

        sub_X = X_shared[union_rows]
        leg_vec = legacy_positive_rank_scores(sub_X, shared_idx, {prog: keep}, max_rank=1500)[prog].to_numpy()
        if requires_debiased:
            c_spec = {prog: controls_map[(lin, prog)]}
            deb_df, raw_u_df, _ = debiased_ucell_scores(
                sub_X, shared_idx, {prog: keep}, matched_controls=c_spec, max_rank=1500, return_components=True
            )
            log_deb_vec = log1p_cp10k_scores(
                sub_X, shared_idx, {prog: keep}, matched_controls=c_spec, input_scale="counts"
            )[prog].to_numpy()
            deb_vec = deb_df[prog].to_numpy()
            raw_u_vec = raw_u_df[prog].to_numpy()
            primary_vec = deb_vec

            filt = [g for g in keep if g not in LIVER_RBC and g not in ("HLA-DRB3", "HLA-DRB4", "HLA-DRB5")]
            if len(filt) >= 5:
                m3_var_key = (lin, f"{prog}__remove_liver_rbc_HLA_DRB345")
                nl_ctrl = controls_map.get(m3_var_key, controls_map[(lin, prog)])
                nl_vec = debiased_ucell_scores(
                    sub_X, shared_idx, {prog: filt}, matched_controls={prog: nl_ctrl}, max_rank=1500
                )[prog].to_numpy()
            else:
                nl_vec = np.full(len(union_rows), np.nan)
        else:
            raw_u_vec = ucell_scores(sub_X, shared_idx, {prog: keep}, max_rank=1500)[prog].to_numpy()
            log_deb_vec = log1p_cp10k_scores(
                sub_X, shared_idx, {prog: keep}, debiased=False, input_scale="counts"
            )[prog].to_numpy()
            deb_vec = leg_vec
            primary_vec = leg_vec
            nl_vec = leg_vec

        cell_score_frames.append(
            pd.DataFrame(
                {
                    "sample_uid": sample_arr[union_rows],
                    "donor_uid": obs_all["donor_uid"].to_numpy(str)[union_rows],
                    "tissue": obs_all["tissue"].to_numpy(str)[union_rows],
                    "treatment_status": obs_all["treatment_status"].to_numpy(str)[union_rows],
                    "lineage": lin,
                    "program": prog,
                    "subtype_m2": sub[union_rows],
                    "myeloid_fine": mf[union_rows],
                    "in_main": m_primary[union_rows],
                    "in_b03": m_b03[union_rows],
                    "in_g500": m_g500[union_rows],
                    "in_g800": m_g800[union_rows],
                    "in_g200": m_g200[union_rows],
                    "in_nomt": m_nomt[union_rows],
                    "score": primary_vec,
                    "ucell_debiased": deb_vec,
                    "ucell_raw": raw_u_vec,
                    "log1p_cp10k_debiased": log_deb_vec,
                    "legacy_positive": leg_vec,
                    "drop_liver_Hb_DRB345": nl_vec,
                }
            )
        )

    scored_cells = (
        pd.concat(cell_score_frames, ignore_index=True)
        if cell_score_frames
        else pd.DataFrame(
            columns=[
                "sample_uid", "donor_uid", "tissue", "treatment_status", "lineage", "program",
                "subtype_m2", "myeloid_fine", "in_main", "in_b03", "in_g500", "in_g800",
                "in_g200", "in_nomt", "score", "ucell_debiased", "ucell_raw",
                "log1p_cp10k_debiased", "legacy_positive", "drop_liver_Hb_DRB345",
            ]
        )
    )

    def build_sample_table(mask_col: str, val_col: str) -> pd.DataFrame:
        sub_df = scored_cells.loc[scored_cells[mask_col]] if not scored_cells.empty else scored_cells
        agg = (
            sub_df.groupby(["sample_uid", "lineage", "program"])[val_col].agg(["mean", "size"]).reset_index()
            if not sub_df.empty
            else pd.DataFrame(columns=["sample_uid", "lineage", "program", "mean", "size"])
        )
        agg_idx = agg.set_index(["sample_uid", "lineage", "program"])
        rows = []
        for r in eval_registry.itertuples(index=False):
            for s in samples.itertuples(index=False):
                key = (s.sample_uid, str(r.lineage), str(r.program))
                if key in agg_idx.index:
                    n_c = int(agg_idx.at[key, "size"])
                    m_v = float(agg_idx.at[key, "mean"]) if n_c > 0 else np.nan
                else:
                    n_c = 0
                    m_v = np.nan
                rows.append(
                    {
                        "sample_uid": s.sample_uid,
                        "lineage": str(r.lineage),
                        "program": str(r.program),
                        "n_cells": n_c,
                        "unfiltered_score": m_v,
                        "score": m_v,
                    }
                )
        return pd.DataFrame(rows)

    policies = [
        ("main", "in_main", "score", 50),
        ("n20", "in_main", "score", 20),
        ("n100", "in_main", "score", 100),
        ("rule_annotation_m1", "in_main", "score", 50),
        ("rule_annotation_m1_n20", "in_main", "score", 20),
        ("rule_annotation_m1_log1p_debiased", "in_main", "log1p_cp10k_debiased", 50),
        ("author_b03_500g_15mt", "in_b03", "score", 50),
        ("author_b03_500g_15mt_n20", "in_b03", "score", 20),
        ("min_genes_500", "in_g500", "score", 50),
        ("min_genes_800_align164", "in_g800", "score", 50),
        ("min_genes_200", "in_g200", "score", 50),
        ("no_mt_filter", "in_nomt", "score", 50),
        ("log1p_cp10k_debiased", "in_main", "log1p_cp10k_debiased", 50),
        ("ucell_raw", "in_main", "ucell_raw", 50),
        ("legacy_positive", "in_main", "legacy_positive", 50),
        ("drop_liver_Hb_DRB345", "in_main", "drop_liver_Hb_DRB345", 50),
    ]
    sens_summaries = []
    policy_outputs: dict[str, dict[str, pd.DataFrame]] = {}
    for pol_name, m_col, v_col, cutoff in policies:
        stab = build_sample_table(m_col, v_col)
        stab.loc[stab["n_cells"] < cutoff, "score"] = np.nan
        pat, les, tis, _ = calculate_patient_contrasts(samples, stab[SCORE_COLUMNS])
        summ = summarize_patients(pat, seed=SEED)
        policy_outputs[pol_name] = {
            "sample_scores": stab,
            "patient_contrasts": pat,
            "effects_summary": summ,
            "patient_tissue_means": tis,
            "lesion_means": les,
        }
        summ_tag = summ.copy()
        summ_tag.insert(0, "analysis", pol_name)
        sens_summaries.append(summ_tag)

    # Focal fine-subtype contrasts (`Macrophage`, `CD4 T cell`, etc.) at `>=50` and `>=20`
    focal_summaries, focal_patients = [], []
    main_cells = scored_cells.loc[scored_cells["in_main"]] if not scored_cells.empty else scored_cells
    for min_cells in (50, 20):
        for r in eval_registry.itertuples(index=False):
            lin, prog = str(r.lineage), str(r.program)
            parent_lin = PARENT_LINEAGE_MAP.get(lin, lin)
            targets = FOCAL_SUBTYPES_BY_PARENT.get(parent_lin, ())
            p_rows = main_cells.loc[main_cells["lineage"].eq(lin) & main_cells["program"].eq(prog)] if not main_cells.empty else main_cells
            for st in targets:
                if p_rows.empty:
                    s_agg = pd.DataFrame(columns=["sample_uid", "mean", "size"])
                else:
                    if st == "monocyte/macrophage":
                        smask = p_rows["subtype_m2"].eq("monocyte/macrophage")
                    elif st in ("Macrophage", "Monocyte"):
                        smask = p_rows["subtype_m2"].eq("monocyte/macrophage") & p_rows["myeloid_fine"].eq(st)
                    else:
                        smask = p_rows["subtype_m2"].eq(st)
                    s_agg = p_rows.loc[smask].groupby("sample_uid")["score"].agg(["mean", "size"]).reset_index()
                s_tab = pd.DataFrame({"sample_uid": samples["sample_uid"], "lineage": st, "program": prog}).merge(
                    s_agg, on="sample_uid", how="left"
                )
                s_tab["n_cells"] = s_tab["size"].fillna(0).astype(int)
                s_tab["score"] = s_tab["mean"].where(s_tab["n_cells"] >= min_cells)
                p_df, _, _, _ = calculate_patient_contrasts(samples, s_tab[SCORE_COLUMNS])
                s_df = summarize_patients(p_df, seed=SEED)
                p_df.insert(0, "min_cells_threshold", min_cells)
                p_df.insert(0, "parent_candidate_lineage", lin)
                p_df.insert(0, "annotation_branch", "rule_annotation_m1")
                s_df.insert(0, "min_cells_threshold", min_cells)
                s_df.insert(0, "parent_candidate_lineage", lin)
                s_df.insert(0, "annotation_branch", "rule_annotation_m1")
                focal_patients.append(p_df)
                focal_summaries.append(s_df)

    # Treatment-naive (COL07, COL12, COL16) vs Preoperative Chemotherapy (COL15, COL17, COL18) descriptive stratification
    donor_tx_map = samples[["donor_uid", "treatment_status", "preop_regimen"]].drop_duplicates().set_index("donor_uid")
    tx_summ_rows = []
    tx_pat_frames = []
    for cutoff_label, pol_key in (("n50_primary", "main"), ("n20_support", "n20")):
        pat_base = policy_outputs[pol_key]["patient_contrasts"].copy()
        pat_base["treatment_status"] = pat_base["donor_uid"].map(donor_tx_map["treatment_status"])
        pat_base["preop_regimen"] = pat_base["donor_uid"].map(donor_tx_map["preop_regimen"])
        pat_base.insert(0, "support_policy", cutoff_label)
        tx_pat_frames.append(pat_base)
        for tx_group in ("treatment_naive", "preoperative_chemotherapy"):
            sub_pat = pat_base.loc[pat_base["treatment_status"].eq(tx_group)]
            if sub_pat.empty:
                continue
            sub_sum = summarize_patients(sub_pat.drop(columns=["support_policy", "treatment_status", "preop_regimen"]), seed=SEED)
            sub_sum = sub_sum.loc[sub_sum["contrast"].eq("D_T")].copy()
            # Attach min/max D_T for small-N transparency
            minmax = sub_pat.groupby(["lineage", "program"])["D_T"].agg(min_D_T="min", max_D_T="max").reset_index()
            sub_sum = sub_sum.merge(minmax, on=["lineage", "program"], how="left")
            sub_sum.insert(0, "treatment_status", tx_group)
            sub_sum.insert(0, "support_policy", cutoff_label)
            sub_sum["inference_scope"] = "descriptive_small_subgroup_N_le_3_no_causal_treatment_claim"
            tx_summ_rows.append(sub_sum)

    # Symmetric subtype composition vs within-subtype state decomposition (incl. B vs Plasma UPR)
    if not main_cells.empty:
        decomp_cells = main_cells.copy()
        decomp_cells["fine_sub"] = np.where(
            decomp_cells["lineage"].eq("Mono/Mac/DC")
            & decomp_cells["subtype_m2"].eq("monocyte/macrophage")
            & decomp_cells["myeloid_fine"].isin(["Macrophage", "Monocyte"]),
            decomp_cells["myeloid_fine"],
            decomp_cells["subtype_m2"],
        )
        tot_per_smp = decomp_cells.groupby(["sample_uid", "lineage", "program"]).size().rename("n_cells").reset_index()
        sub_per_smp = (
            decomp_cells.groupby(["sample_uid", "lineage", "program", "fine_sub"])["score"]
            .agg(mean_score="mean", n_subtype="size")
            .reset_index()
            .rename(columns={"fine_sub": "subtype"})
        )
        subtype_sample_df = sub_per_smp.merge(tot_per_smp, on=["sample_uid", "lineage", "program"], how="left")
        subtype_decomp = decompose_subtypes(samples, subtype_sample_df)
        bp_decomp = subtype_decomp.loc[subtype_decomp["lineage"].eq("B/Plasma")].reset_index(drop=True)
    else:
        subtype_decomp = pd.DataFrame()
        bp_decomp = pd.DataFrame()

    return {
        "sample_immune_composition": pd.DataFrame(comp_rows),
        "sample_lineage_support": pd.DataFrame(support_rows),
        "candidate_measurement": pd.DataFrame(measured_rows),
        "sensitivity_summary": pd.concat(sens_summaries, ignore_index=True),
        "focal_subtype_effects_summary": pd.concat(focal_summaries, ignore_index=True) if focal_summaries else pd.DataFrame(),
        "focal_subtype_patient_contrasts": pd.concat(focal_patients, ignore_index=True) if focal_patients else pd.DataFrame(),
        "treatment_stratified_summary": pd.concat(tx_summ_rows, ignore_index=True) if tx_summ_rows else pd.DataFrame(),
        "treatment_stratified_patient_contrasts": pd.concat(tx_pat_frames, ignore_index=True) if tx_pat_frames else pd.DataFrame(),
        "patient_subtype_decomposition": subtype_decomp,
        "patient_B_plasma_decomposition": bp_decomp,
        "policy_outputs": policy_outputs,
    }


def run(
    freeze_dir: Path = FROZEN,
    out_dir: Path = OUT,
    processed_dir: Path = PROCESSED,
) -> dict[str, object]:
    t0 = time.perf_counter()
    started_utc = utc_now()
    freeze_dir = Path(freeze_dir).resolve()
    out_dir = Path(out_dir).resolve()
    processed_dir = Path(processed_dir).resolve()
    root_resolved = ROOT.resolve()
    _ = processed_dir.relative_to(root_resolved)
    if processed_dir in (
        (root_resolved / "data/processed/m1").resolve(),
        (root_resolved / "data/processed/m2").resolve(),
    ):
        raise ValueError(f"Refusing to write external GSE178318 H5ADs into development directory: {processed_dir}")

    verify_freeze_integrity(freeze_dir)
    registry = pd.read_csv(freeze_dir / "frozen_candidate_registry.csv")
    if registry.empty or not registry["status"].astype(str).str.lower().eq("frozen").all():
        raise ValueError("A frozen M5 candidate registry is required")

    sets = pd.read_csv(freeze_dir / "gene_sets/program_genes.csv")
    model_genes_path = freeze_dir / "gene_sets/model_program_genes.csv"
    sec_path = freeze_dir / "secondary_sensitivity_candidates.csv"
    eval_registry = registry.copy()
    if model_genes_path.is_file():
        mg_df = pd.read_csv(model_genes_path)
        sets = pd.concat(
            [sets[["lineage", "program", "gene_symbol"]], mg_df[["lineage", "program", "gene_symbol"]]],
            ignore_index=True,
        ).drop_duplicates(["lineage", "program", "gene_symbol"])
        extra_pairs: list[tuple[str, str]] = []
        if sec_path.is_file():
            sec_df = pd.read_csv(sec_path)
            extra_pairs.extend(zip(sec_df["lineage"].astype(str), sec_df["program"].astype(str)))
        extra_pairs.append(("B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"))
        mg_keys = set(zip(mg_df["lineage"].astype(str), mg_df["program"].astype(str)))
        existing = set(zip(eval_registry["lineage"].astype(str), eval_registry["program"].astype(str)))
        add_rows = []
        template = eval_registry.iloc[0].to_dict()
        for lin, prog in extra_pairs:
            if (lin, prog) in mg_keys and (lin, prog) not in existing:
                existing.add((lin, prog))
                r_copy = dict(template)
                r_copy.update(
                    candidate_id=f"{lin}::{prog}",
                    lineage=lin,
                    program=prog,
                    min_measured_positive_genes=5,
                    min_measured_fraction_positive=0.8,
                )
                add_rows.append(r_copy)
        if add_rows:
            eval_registry = pd.concat([eval_registry, pd.DataFrame(add_rows)], ignore_index=True)

    shared_axis = list(load_shared_nuclear_universe(SHARED_UNIVERSE_NAME))
    cand_keys = sorted(set(zip(eval_registry["lineage"].astype(str), eval_registry["program"].astype(str))))
    variant_keys = [(lin, f"{prog}__remove_liver_rbc_HLA_DRB345") for lin, prog in cand_keys]
    controls_map = load_frozen_matched_controls(
        freeze_dir, cand_keys + variant_keys, registry, shared_axis=shared_axis, strict_keys=cand_keys
    )

    barcodes_path = RAW / "GSE178318_barcodes.tsv.gz"
    genes_path = RAW / "GSE178318_genes.tsv.gz"
    matrix_path = RAW / "GSE178318_matrix.mtx.gz"
    for req_file in (barcodes_path, genes_path, matrix_path):
        if not req_file.is_file():
            raise FileNotFoundError(f"Missing GSE178318 raw input file: {req_file}")

    with gzip.open(barcodes_path, "rt", encoding="utf-8") as handle:
        barcodes = [line.rstrip("\r\n") for line in handle if line.strip()]
    cell_df, gsm_audit_df = parse_and_audit_barcodes(barcodes)

    genes_raw = pd.read_csv(
        genes_path,
        sep="\t",
        header=None,
        names=["source_feature_id", "gene_symbol"],
        dtype=str,
        keep_default_na=False,
    )
    if genes_raw["source_feature_id"].duplicated().any() or (genes_raw == "").any().any():
        raise ValueError("Invalid or duplicate Ensembl IDs in GSE178318_genes.tsv.gz")

    excl_353 = set(load_shared_nuclear_universe("excluded_353_canonical"))
    shared_164_set = set(load_shared_nuclear_universe("245_164_shared_nuclear"))
    shared_203_set = set(shared_axis)
    norm_genes = normalize_gene_axis(
        genes_raw["gene_symbol"].tolist(),
        genes_raw["source_feature_id"].tolist(),
        excluded_canonical=excl_353,
    )
    raw_mt_mask = norm_genes["canonical_symbol"].astype(str).str.startswith("MT-").to_numpy(dtype=bool)
    measurable_mask = norm_genes["is_unique_measurable"].to_numpy(dtype=bool)

    var_measurable = pd.DataFrame(
        {
            "source_feature_id": genes_raw.loc[measurable_mask, "source_feature_id"].to_numpy(dtype=str),
            "gene_symbol": genes_raw.loc[measurable_mask, "gene_symbol"].to_numpy(dtype=str),
            "canonical_symbol": norm_genes.loc[measurable_mask, "canonical_symbol"].to_numpy(dtype=str),
            "hgnc_id": norm_genes.loc[measurable_mask, "hgnc_id"].to_numpy(dtype=str),
            "hgnc_match_type": norm_genes.loc[measurable_mask, "match_type"].to_numpy(dtype=str),
            "in_excluded_353": norm_genes.loc[measurable_mask, "in_excluded_353"].to_numpy(dtype=bool),
        }
    )
    var_measurable["in_shared_nuclear_245_164"] = var_measurable["canonical_symbol"].isin(shared_164_set).to_numpy()
    var_measurable["in_shared_nuclear_245_164_203608"] = var_measurable["canonical_symbol"].isin(shared_203_set).to_numpy()
    var_measurable.index = pd.Index(var_measurable["source_feature_id"], name="feature_id")

    al_shared = align_axis_to_canonical_targets(
        genes_raw["gene_symbol"].tolist(), shared_axis, genes_raw["source_feature_id"].tolist()
    )
    if not al_shared["present"].all():
        raise ValueError("GSE178318 missing expected genes from 245_164_203608_shared_nuclear")

    # Map only the 12 tumor PT/LM GSMs into memory (hard physical exclusion of the 3 PBMC GSMs)
    tumor_mask = cell_df["is_tumor_pt_lm"].to_numpy(dtype=bool)
    tumor_cell_map = np.full(len(cell_df), -1, dtype=np.int32)
    tumor_cell_map[tumor_mask] = np.arange(int(tumor_mask.sum()), dtype=np.int32)
    tumor_cells_df = cell_df.loc[tumor_mask].reset_index(drop=True)

    X_tumor_raw, nnz_all, total_all_umi = stream_tumor_count_matrix(
        matrix_path, len(genes_raw), len(cell_df), tumor_cell_map
    )

    processed_dir.mkdir(parents=True, exist_ok=True)
    for pbmc_gsm in FORBIDDEN_PBMC_GSMS:
        pbmc_file = processed_dir / f"{pbmc_gsm}.h5ad"
        if pbmc_file.is_file():
            pbmc_file.unlink()

    qc_rows = []
    obs_parts = []
    X_shared_parts = []
    meas_canon_idx = pd.Index(var_measurable["canonical_symbol"].astype(str))
    shared_col_in_meas = meas_canon_idx.get_indexer(shared_axis)
    if (shared_col_in_meas < 0).any():
        raise ValueError("Measurable canonical axis missing genes from 245_164_203608_shared_nuclear")

    tumor_manifest_rows = []
    for meta_row in gsm_audit_df.loc[gsm_audit_df["is_tumor_pt_lm"]].to_dict("records"):
        gsm = str(meta_row["sample_uid"])
        local_idx = np.flatnonzero(tumor_cells_df["sample_uid"].eq(gsm).to_numpy())
        X_sub_raw = X_tumor_raw[local_idx]
        sub_cells = tumor_cells_df.iloc[local_idx].reset_index(drop=True)

        adata, qc_info = annotate_and_qc_sample(
            X_sub_raw, meta_row, sub_cells, raw_mt_mask, measurable_mask, var_measurable
        )
        h5ad_path = processed_dir / f"{gsm}.h5ad"
        rel_h5ad = str(h5ad_path.relative_to(root_resolved))
        reuse_existing = False
        if h5ad_path.is_file():
            existing_obj = ad.read_h5ad(h5ad_path)
            if (
                sparse.isspmatrix_csr(existing_obj.X)
                and existing_obj.shape == adata.shape
                and existing_obj.X.nnz == adata.X.nnz
                and int(existing_obj.X.sum(dtype=np.int64)) == int(adata.X.sum(dtype=np.int64))
                and existing_obj.obs_names.equals(adata.obs_names)
                and existing_obj.var_names.equals(adata.var_names)
                and existing_obj.obs["qc_retained"].astype(bool).equals(adata.obs["qc_retained"].astype(bool))
                and existing_obj.obs["lineage_m2"].astype(str).equals(adata.obs["lineage_m2"].astype(str))
            ):
                reuse_existing = True
                stored_bytes = h5ad_path.stat().st_size
        if not reuse_existing:
            stored_bytes = write_checked_h5ad(adata, h5ad_path)

        # Verify zero-friction compatibility with model_bags.sample_bags & specimen_immune_composition for Baton 2
        _ = sample_bags(h5ad_path, min_cells=50, max_cells=None, gene_universe=SHARED_UNIVERSE_NAME)
        _ = specimen_immune_composition(h5ad_path)

        qc_info["h5ad_path"] = rel_h5ad
        qc_info["h5ad_bytes"] = stored_bytes
        qc_info["h5ad_sha256"] = sha(h5ad_path)
        qc_rows.append(qc_info)

        obs_parts.append(adata.obs.copy())
        X_shared_parts.append(adata.X[:, shared_col_in_meas].tocsr())
        tumor_manifest_rows.append(
            {
                "dataset": DATASET,
                "sample_uid": gsm,
                "donor_uid": str(meta_row["donor_uid"]),
                "lesion_uid": str(meta_row["lesion_uid"]),
                "patient_id": str(meta_row["patient_id"]),
                "barcode_tag": str(meta_row["barcode_tag"]),
                "tissue": str(meta_row["tissue"]),
                "role": "external",
                "treatment_status": str(meta_row["treatment_status"]),
                "preop_regimen": str(meta_row["preop_regimen"]),
                "sex": str(meta_row["sex"]),
                "age": int(meta_row["age"]),
                "stage": str(meta_row["stage"]),
                "msi_status": str(meta_row["msi_status"]),
                "input_scale": "raw_counts",
                "eligible_for_paired_analysis": True,
                "patient_mapping_verified": True,
                "counts_values_verified": True,
                "input_cells": qc_info["input_cells"],
                "retained_cells": qc_info["retained_cells_m2_primary"],
                "h5ad_path": qc_info["h5ad_path"],
                "h5ad_sha256": qc_info["h5ad_sha256"],
            }
        )
        print(
            f"[validate_gse178318] {gsm} ({meta_row['barcode_tag']}, {meta_row['treatment_status']}): "
            f"{qc_info['input_cells']:,} -> {qc_info['retained_cells_m2_primary']:,} retained",
            flush=True,
        )

    samples_df = pd.DataFrame(tumor_manifest_rows)
    obs_all = pd.concat(obs_parts, axis=0)
    X_shared_all = sparse.vstack(X_shared_parts, format="csr")

    res = evaluate_tumor_cohort(
        X_shared_all,
        shared_axis,
        obs_all,
        samples_df,
        eval_registry,
        sets,
        controls_map,
        requires_debiased=registry_requires_debiased(registry),
    )

    coverage_dir = out_dir / "coverage"
    pe_dir = out_dir / "patient_effects"
    reg_dir = out_dir / "registered_tests"
    for d in (out_dir, coverage_dir, pe_dir, reg_dir):
        d.mkdir(parents=True, exist_ok=True)

    gsm_audit_df.to_csv(coverage_dir / "all_15_gsm_barcode_audit.csv", index=False)
    samples_df.to_csv(coverage_dir / "sample_manifest.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(coverage_dir / "qc_and_exclusions.csv", index=False)
    res["sample_immune_composition"].to_csv(coverage_dir / "sample_immune_composition.csv", index=False)
    res["sample_lineage_support"].to_csv(coverage_dir / "sample_lineage_support.csv", index=False)
    res["candidate_measurement"].to_csv(coverage_dir / "candidate_measurement.csv", index=False)
    (coverage_dir / "external_rank_background.txt").write_text("\n".join(shared_axis) + "\n", encoding="utf-8")

    main_out = res["policy_outputs"]["main"]
    main_out["sample_scores"].to_csv(pe_dir / "sample_program_scores.csv", index=False)
    main_out["patient_contrasts"].to_csv(pe_dir / "patient_contrasts.csv", index=False)
    main_out["effects_summary"].to_csv(pe_dir / "effects_summary.csv", index=False)
    main_out["patient_tissue_means"].to_csv(pe_dir / "patient_tissue_means.csv", index=False)
    main_out["lesion_means"].to_csv(pe_dir / "lesion_means.csv", index=False)
    leave_one_out(main_out["patient_contrasts"]).to_csv(pe_dir / "leave_one_patient_out.csv", index=False)
    res["sensitivity_summary"].to_csv(pe_dir / "sensitivity_summary.csv", index=False)
    res["focal_subtype_effects_summary"].to_csv(pe_dir / "focal_subtype_effects_summary.csv", index=False)
    res["focal_subtype_patient_contrasts"].to_csv(pe_dir / "focal_subtype_patient_contrasts.csv", index=False)
    res["treatment_stratified_summary"].to_csv(pe_dir / "treatment_stratified_summary.csv", index=False)
    res["treatment_stratified_patient_contrasts"].to_csv(pe_dir / "treatment_stratified_patient_contrasts.csv", index=False)
    res["patient_subtype_decomposition"].to_csv(pe_dir / "patient_subtype_decomposition.csv", index=False)
    res["patient_B_plasma_decomposition"].to_csv(pe_dir / "patient_B_plasma_decomposition.csv", index=False)

    for pol_save in (
        "n20",
        "n100",
        "author_b03_500g_15mt",
        "min_genes_800_align164",
        "log1p_cp10k_debiased",
        "ucell_raw",
        "legacy_positive",
        "drop_liver_Hb_DRB345",
    ):
        res["policy_outputs"][pol_save]["patient_contrasts"].to_csv(
            pe_dir / f"patients_{pol_save}.csv", index=False
        )

    # Cohort-local Holm adjustment for SECONDARY_DT_178_V1_1
    plan_full = pd.read_csv(freeze_dir / "registered_tests.csv", keep_default_na=False)
    plan_178 = plan_full.loc[plan_full["cohort"].eq(DATASET)].reset_index(drop=True)
    adj_178, expl_178, fam_178, audit_178 = adjust_registered(plan_178, main_out["effects_summary"])
    adj_178.to_csv(reg_dir / "registered_results.csv", index=False)
    fam_178.to_csv(reg_dir / "family_summary.csv", index=False)
    expl_178.to_csv(reg_dir / "unplanned_results.csv", index=False)
    (reg_dir / "audit.json").write_text(json.dumps(audit_178, indent=2) + "\n", encoding="utf-8")

    profile = {
        "dataset": DATASET,
        "validation_role": "authorized_prospective_zero_exposure_secondary_external_cohort",
        "freeze_version": "CRLM_V1.1_245_F1",
        "started_utc": started_utc,
        "completed_utc": utc_now(),
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "barcodes_sha256": sha(barcodes_path),
        "genes_sha256": sha(genes_path),
        "matrix_sha256": sha(matrix_path),
        "code_sha256": sha(Path(__file__)),
        "program_scores_sha256": sha(ROOT / "src/crlm/program_scores.py"),
        "annotation_m1_sha256": sha(ROOT / "src/crlm/annotation_m1.py"),
        "hgnc_alias_map_sha256": sha(ROOT / "src/crlm/hgnc_alias_map.py"),
        "all_15_gsm_cells": int(len(cell_df)),
        "excluded_pbmc_gsms": 3,
        "excluded_pbmc_cells": int((~tumor_mask).sum()),
        "tumor_pt_lm_gsms": int(len(samples_df)),
        "paired_donors": int(samples_df["donor_uid"].nunique()),
        "treatment_naive_donors": 3,
        "preoperative_chemotherapy_donors": 3,
        "tumor_input_cells": int(len(obs_all)),
        "tumor_retained_cells_m2_primary": int(obs_all["qc_retained"].sum()),
        "tumor_retained_cells_author_b03_500g_15mt": int(obs_all["qc_author_b03_500g_15mt"].sum()),
        "source_nnz_all_15": int(nnz_all),
        "source_total_umi_all_15": int(total_all_umi),
        "raw_genes": int(len(genes_raw)),
        "unique_measurable_canonical_genes": int(len(var_measurable)),
        "shared_nuclear_axis_name": SHARED_UNIVERSE_NAME,
        "shared_nuclear_axis_genes": int(len(shared_axis)),
    }
    (out_dir / "run_profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
    print(adj_178[["test_id", "lineage", "program", "n_patients", "n_positive", "n_negative", "mean_difference", "raw_p", "holm_adjusted_p", "direction_status"]].to_string(index=False), flush=True)
    return profile


def run_selfcheck() -> None:
    """Pure in-memory & TemporaryDirectory verification (zero project disk writes)."""
    rng = np.random.default_rng(SEED)
    full_shared = list(load_shared_nuclear_universe(SHARED_UNIVERSE_NAME))
    full_shared_set = set(full_shared)

    # 1. Verify barcode audit on real downloaded barcodes if present
    bpath = RAW / "GSE178318_barcodes.tsv.gz"
    if bpath.is_file():
        with gzip.open(bpath, "rt", encoding="utf-8") as fh:
            bcs = [ln.strip() for ln in fh if ln.strip()]
        cdf, gdf = parse_and_audit_barcodes(bcs)
        assert len(cdf) == 140281 and len(gdf) == 15
        assert int(gdf.loc[gdf["is_tumor_pt_lm"], "observed_cells"].sum()) == 113331
        assert int(gdf.loc[~gdf["is_tumor_pt_lm"], "observed_cells"].sum()) == 26950

    # 2. Verify freeze_dir resolution & missing-control guard
    with tempfile.TemporaryDirectory() as tmp:
        tdir = Path(tmp)
        (tdir / "program_scores.py").write_bytes((ROOT / "src/crlm/program_scores.py").read_bytes())
        (tdir / "gene_sets").mkdir()
        mc_file = tdir / "gene_sets/matched_controls.csv"
        pd.DataFrame(
            [{"lineage": "Mono/Mac/DC", "program": "SYN_SAM", "target_gene": "APOE", "replicate": 0, "control_gene": "A1BG"}]
        ).to_csv(mc_file, index=False)
        (tdir / "frozen_inputs_sha256.txt").write_text(
            f"{sha(mc_file)}  artifacts/m5/gene_sets/matched_controls.csv\n", encoding="utf-8"
        )
        assert verify_freeze_integrity(Path(os.path.relpath(tdir, Path.cwd()))) == 1
        reg_bad = pd.DataFrame(
            [
                {"lineage": "Mono/Mac/DC", "program": "SYN_SAM", "scoring_method": "ucell_debiased"},
                {"lineage": "T/NK", "program": "MISSING", "scoring_method": "ucell_debiased"},
            ]
        )
        raised = False
        try:
            load_frozen_matched_controls(tdir, [("Mono/Mac/DC", "SYN_SAM"), ("T/NK", "MISSING")], reg_bad, shared_axis=full_shared)
        except ValueError:
            raised = True
        assert raised

        # 3. Synthetic end-to-end M2 QC + annotation_m1 + H5AD model_bags check + treatment stratification
        syn_axis = full_shared[:450]
        inject = [
            "CD3D", "CD3E", "CD4", "CD8A", "TRAC", "IL7R",
            "MS4A1", "CD79A", "MZB1", "JCHAIN", "SDC1", "TNFRSF17",
            "C1QA", "C1QB", "CD14", "CD68", "CSF1R", "SPP1", "PTPRC", "EPCAM", "KRT19",
        ]
        for idx_m, g_m in enumerate(inject):
            if g_m not in syn_axis:
                syn_axis[idx_m] = g_m
        raw_axis = syn_axis + ["MT-CO1", "MT-ND1"]
        var_syn = pd.DataFrame(
            {
                "source_feature_id": [f"ENSG{i:011d}" for i in range(len(raw_axis))],
                "gene_symbol": raw_axis,
                "canonical_symbol": raw_axis,
                "hgnc_id": [f"HGNC:{i+1}" for i in range(len(raw_axis))],
                "hgnc_match_type": "approved_symbol",
                "in_excluded_353": [g.startswith("MT-") for g in raw_axis],
                "in_shared_nuclear_245_164": [not g.startswith("MT-") for g in raw_axis],
                "in_shared_nuclear_245_164_203608": [not g.startswith("MT-") for g in raw_axis],
            }
        )
        var_syn.index = pd.Index(var_syn["source_feature_id"], name="feature_id")
        all_panel = {g for _, a, s in PANELS.values() for g in (a + " " + s).split()} | set(AUXILIARY_MARKERS)

        scale_row = rng.uniform(1.5, 6.0, size=(122, 1))
        X_mat = rng.poisson(scale_row, size=(122, len(raw_axis))).astype(np.int32)
        for col_i, g_name in enumerate(raw_axis):
            if g_name in all_panel:
                X_mat[:, col_i] = 0
        g_idx = {g: i for i, g in enumerate(raw_axis)}
        for gm in ("PTPRC", "C1QA", "C1QB", "CD68", "CSF1R"):
            X_mat[:70, g_idx[gm]] = 18
        for gm in ("PTPRC", "CD3D", "CD3E", "CD4", "TRAC", "IL7R"):
            X_mat[70:120, g_idx[gm]] = 18
        X_mat[:120, g_idx["MT-CO1"]] = 5
        # Cell 120: high MT (>20%); Cell 121: below 300-gene floor
        X_mat[120, :] = X_mat[0, :]
        X_mat[120, g_idx["MT-CO1"]] = 1500
        X_mat[121, :] = 0
        X_mat[121, :50] = 2

        scells = pd.DataFrame(
            {
                "deposited_barcode": [f"BC{i:04d}_COL07_CRC" for i in range(122)],
                "barcode_16bp": [f"BC{i:04d}" for i in range(122)],
            }
        )
        meta = {
            "sample_uid": "GSM5387660",
            "donor_uid": "GSE178318:COL07",
            "lesion_uid": "GSE178318:COL07_CRC",
            "patient_id": "COL07",
            "tissue": "PT",
            "treatment_status": "treatment_naive",
            "preop_regimen": "None",
        }
        adata_syn, qc_syn = annotate_and_qc_sample(
            sparse.csr_matrix(X_mat),
            meta,
            scells,
            np.array([g.startswith("MT-") for g in raw_axis]),
            np.ones(len(raw_axis), dtype=bool),
            var_syn,
        )
        assert qc_syn["retained_cells_m2_primary"] == 120
        assert qc_syn["excluded_high_pct_mt_gt20"] == 1
        assert qc_syn["excluded_below_min_300g_500c"] == 1
        syn_h5ad = tdir / "GSM5387660.h5ad"
        write_checked_h5ad(adata_syn, syn_h5ad)
        bags = sample_bags(syn_h5ad, min_cells=50, max_cells=None, gene_universe=syn_axis)
        assert "Mono/Mac/DC" in bags and "T/NK" in bags and bags["Mono/Mac/DC"].population == 70

    print("validate_gse178318.py --selfcheck PASSED (PBMC isolation, M2 QC, annotation_m1, H5AD bag readback verified).")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selfcheck", action="store_true", help="Run in-memory + temporary self-check.")
    parser.add_argument("--freeze-dir", type=Path, default=FROZEN)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument("--processed-dir", type=Path, default=PROCESSED)
    args = parser.parse_args()
    if args.selfcheck:
        run_selfcheck()
    else:
        run(freeze_dir=args.freeze_dir, out_dir=args.out_dir, processed_dir=args.processed_dir)


if __name__ == "__main__":
    main()
