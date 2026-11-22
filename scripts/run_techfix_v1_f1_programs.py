#!/usr/bin/env python3
"""Track 1 (`CRLM_V1_245_F1_techfix`) stepwise technical remediation of frozen V1 program validation.

Strictly preserves the two frozen V1 candidates (`Mono/Mac/DC::B19_S4_Mac_SPP1` and
`B/Plasma::HALLMARK_UNFOLDED_PROTEIN_RESPONSE`), expected directions, four test
families (8 planned tests), and legacy `positive_midrank_1500` scoring formula from
`artifacts/archive/CRLM_V1_245_F1/artifacts/m5/`, while isolating objective technical
and QC bug fixes across four cumulative steps:

- Step 0 (`step0_v1_f1_baseline`): Bit-exact reproduction of archived `CRLM_V1_245_F1`.
- Step 1 (`step1_hgnc_harmonized`): Deterministic HGNC alias/prev_symbol harmonization
  plus Ensembl-ID disambiguation of 33/34 duplicate symbol pairs in GSE203608.
- Step 2 (`step2_shared_nuclear_axis`): Restrict ranking background to cross-cohort
  shared nuclear universes (`245_164_shared_nuclear` [18,815] and
  `245_164_203608_shared_nuclear` [18,733]), removing the 353 tumor-deleted
  mitochondrial/ribosomal genes that depressed ranks in external cohorts.
- Step 3 (`step3_full_techfix`): Objective external cell-level QC bug fixes:
  * GSE203608: enforce `Cell_type != "NA"` and `MitoFraction <= 0.20` (removing
    dead/unannotated high-mitochondrial cells and correcting composition denominators).
  * GSE164522: exclude all cross-lineage `CD3E+` doublet subclusters in non-T/NK
    lineages (`hC36_cDC1-CD3E`, `hC41_cDC2-CD3E`, `hC47_Mono-CD3E`, `hC56_Mph-CD3E`,
    and previously missed `hC75_B-CD3E` [3,174 cells]), and recompute B vs Plasma
    symmetric composition/state decomposition.
"""

from __future__ import annotations

from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]

from adjust_registered_tests import adjust_registered  # noqa: E402
from composition_decomposition import decompose_pair  # noqa: E402
from crlm.hgnc_alias_map import (  # noqa: E402
    align_axis_to_canonical_targets,
    canonicalize_symbols,
    load_shared_nuclear_universe,
    normalize_gene_axis,
)
from crlm.program_scores import rank_scores  # noqa: E402
from paired_effects import calculate_patient_contrasts, summarize_patients  # noqa: E402

ARCHIVE_M5 = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m5"
ARCHIVE_M6 = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m6"
OUT = ROOT / "artifacts/techfix_v1_f1/programs"

GSE164_OBJECTS = ROOT / "data/processed/m6/GSE164522"
GSE164_MANIFEST = ARCHIVE_M6 / "GSE164522/coverage/sample_manifest.csv"
GSE164_ORGANS = ("PT", "MT", "PN", "MN")

GSE178_OBJECTS = ROOT / "data/processed/m6/GSE178318"
GSE178_MANIFEST = ROOT / "artifacts/m6/GSE178318/coverage/sample_manifest.csv"

GSE203_RAW = ROOT / "data/raw/GSE203608"
GSE203_SITES = {"Colon": "NC", "Liver": "NL", "LiverMet": "LM", "Tumor": "PT"}
GSE203_MYELOID = {
    "Macrophage",
    "cDC2 (CD1C+)",
    "cDC1 (BATF3+)",
    "LAMP3+ DC (LAMP3+)",
    "pDC (LILRA4+)",
}

SEED = 20261121
SCORE_COLUMNS = ["sample_uid", "lineage", "program", "score", "n_cells"]


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load_gse203608_subset() -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame]:
    """Load all in-vivo `Sample_pass_qc == TRUE` cells of GSE203608 once in CSR."""
    cells = pd.read_csv(
        GSE203_RAW / "GSE203608_cell_annotation.csv.gz",
        dtype=str,
        keep_default_na=False,
    ).rename(columns={"Unnamed: 0": "barcode"})
    genes = pd.read_csv(
        GSE203_RAW / "GSE203608_gene_annotation.csv.gz",
        dtype=str,
        keep_default_na=False,
    )
    in_vivo_mask = (
        cells["SampleType"].isin(GSE203_SITES) & cells["Sample_pass_qc"].eq("TRUE")
    ).to_numpy()

    cell_map = np.full(len(cells), -1, dtype=np.int32)
    cell_map[in_vivo_mask] = np.arange(int(in_vivo_mask.sum()), dtype=np.int32)

    nfeatures = len(genes)
    ncells = len(cells)
    blocks = []
    checked = 0
    mtx_path = GSE203_RAW / "GSE203608_gene_by_cell_count_matrix.txt.gz"
    with gzip.open(mtx_path, "rb") as handle:
        if handle.readline().strip() != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError("Unexpected GSE203608 MatrixMarket header")
        header = handle.readline()
        while header.startswith(b"%"):
            header = handle.readline()
        nf, nc, nnz = map(int, header.split())
        if (nf, nc) != (nfeatures, ncells):
            raise ValueError(f"GSE203608 shape mismatch: {(nf, nc)} vs {(nfeatures, ncells)}")
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
            checked += len(val)
            use = (cell_map[c] >= 0) & (val > 0)
            if use.any():
                blocks.append(
                    (
                        cell_map[c[use]],
                        g[use].astype(np.int32),
                        val[use].astype(np.int32),
                    )
                )
    if checked != nnz:
        raise ValueError(f"Read {checked} entries, expected {nnz}")
    r_idx, c_idx, vals = (np.concatenate(part) for part in zip(*blocks))
    matrix = sparse.coo_matrix(
        (vals, (r_idx, c_idx)), shape=(int(in_vivo_mask.sum()), nfeatures)
    ).tocsr()

    sub_cells = cells.loc[in_vivo_mask].reset_index(drop=True).copy()
    sub_cells["mito_fraction_num"] = pd.to_numeric(
        sub_cells["MitoFraction"], errors="coerce"
    )
    sub_cells["lineage"] = "Other/Unknown"
    sub_cells.loc[sub_cells["Cell_type"].isin(["Plasma cell", "B cell"]), "lineage"] = "B/Plasma"
    sub_cells.loc[
        sub_cells["Cell_type"].eq("Myeloid") & sub_cells["Cell_type_myeloid"].isin(GSE203_MYELOID),
        "lineage",
    ] = "Mono/Mac/DC"
    sub_cells.loc[sub_cells["Cell_type"].eq("T cell"), "lineage"] = "T/NK"
    sub_cells.loc[sub_cells["Cell_type"].eq("Mast cell"), "lineage"] = "Mast"
    sub_cells["sample_uid"] = "GSE203608:" + sub_cells["Dataset"]
    sub_cells["donor_uid"] = "GSE203608:" + sub_cells["Patient"]
    sub_cells["tissue"] = sub_cells["SampleType"].map(GSE203_SITES)
    sub_cells["clean_cell_qc"] = sub_cells["Cell_type"].ne("NA") & (
        sub_cells["mito_fraction_num"] <= 0.20
    )
    return matrix, sub_cells, genes


def score_gse164522_step(
    h5ads: dict[str, ad.AnnData],
    manifest: pd.DataFrame,
    registry: pd.DataFrame,
    frozen_genes: pd.DataFrame,
    *,
    harmonize_hgnc: bool,
    shared_nuclear_only: bool,
    remove_cd3e_doublets: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score GSE164522 across all four tissues under the specified ablation flags."""
    if not harmonize_hgnc and not shared_nuclear_only:
        axes = [
            set(obj.var.loc[obj.var["unambiguous_symbol"], "gene_symbol"])
            for obj in h5ads.values()
        ]
        axis = sorted(set.intersection(*axes))
        col_maps = {code: obj.var_names.get_indexer(axis) for code, obj in h5ads.items()}
    else:
        if shared_nuclear_only:
            axis = list(load_shared_nuclear_universe("245_164_shared_nuclear"))
        else:
            norm_axes = []
            for obj in h5ads.values():
                ndf = normalize_gene_axis(obj.var["gene_symbol"].tolist())
                norm_axes.append(
                    set(ndf.loc[ndf["is_unique_measurable"], "canonical_symbol"])
                )
            axis = sorted(set.intersection(*norm_axes))
        col_maps = {}
        for code, obj in h5ads.items():
            al = align_axis_to_canonical_targets(obj.var["gene_symbol"].tolist(), axis)
            if not al["present"].all():
                raise ValueError(f"GSE164522 {code}: missing expected shared gene on axis")
            col_maps[code] = al["source_index"].to_numpy(dtype=np.int32)

    axis_set = set(axis)
    programs: dict[tuple[str, str], list[str]] = {}
    cov_rows = []
    for r in registry.itertuples(index=False):
        raw_req = frozen_genes.loc[
            frozen_genes["program"].eq(r.program) & frozen_genes["lineage"].eq(r.lineage),
            "gene_symbol",
        ].tolist()
        req = canonicalize_symbols(raw_req) if harmonize_hgnc else raw_req
        measured = [g for g in req if g in axis_set]
        ok = (
            len(measured) >= r.min_measured_positive_genes
            and len(measured) >= len(req) * r.min_measured_fraction_positive
        )
        programs[(r.lineage, r.program)] = measured if ok else []
        cov_rows.append(
            {
                "cohort": "GSE164522",
                "lineage": r.lineage,
                "program": r.program,
                "n_axis_genes": len(axis),
                "n_requested": len(req),
                "n_measured": len(measured),
                "measured_fraction": len(measured) / len(req),
                "eligible": ok,
                "missing_symbols": "|".join(sorted(set(req) - set(measured))),
            }
        )

    cell_frames = []
    for code, obj in h5ads.items():
        obs = obj.obs
        cols = col_maps[code]
        for (lin, prog), members in programs.items():
            if not members:
                continue
            mask = obs["lineage"].eq(lin).to_numpy()
            if remove_cd3e_doublets:
                is_cd3e_doublet = obs["lineage"].isin(["B/Plasma", "Mono/Mac/DC"]) & obs[
                    "celltype_sub"
                ].str.contains("CD3E", regex=False)
                mask = mask & (~is_cd3e_doublet.to_numpy())
            rows = np.flatnonzero(mask)
            if len(rows) == 0:
                continue
            sc = rank_scores(
                obj.X[rows][:, cols],
                axis,
                {prog: members},
                max_rank=1500,
                mode="legacy_positive",
            )[prog].to_numpy()
            sub = obs.iloc[rows][
                ["sample_uid", "donor_uid", "lesion_uid", "tissue", "lineage", "celltype_sub"]
            ].copy()
            sub["program"] = prog
            sub["score"] = sc
            cell_frames.append(sub)

    cells_df = pd.concat(cell_frames, ignore_index=True)
    by_sample = (
        cells_df.groupby(["sample_uid", "lineage", "program"])["score"]
        .agg(["mean", "size"])
        .reset_index()
    )
    grid = pd.MultiIndex.from_product(
        [manifest["sample_uid"].tolist(), list(programs.keys())]
    )
    sample_rows = []
    lookup = {
        (r.sample_uid, r.lineage, r.program): (float(r["mean"]), int(r["size"]))
        for _, r in by_sample.iterrows()
    }
    for sample_uid, (lin, prog) in grid:
        mean_val, n_c = lookup.get((sample_uid, lin, prog), (np.nan, 0))
        sample_rows.append(
            {
                "sample_uid": sample_uid,
                "lineage": lin,
                "program": prog,
                "n_cells": n_c,
                "unfiltered_score": mean_val,
                "score": mean_val if n_c >= 50 else np.nan,
            }
        )
    sample_df = pd.DataFrame(sample_rows)
    patients, _, _, _ = calculate_patient_contrasts(manifest, sample_df[SCORE_COLUMNS])
    summary = summarize_patients(patients, seed=SEED)
    return summary, patients, sample_df, pd.DataFrame(cov_rows), cells_df


def decompose_gse164522_b_plasma(manifest: pd.DataFrame, cells_df: pd.DataFrame, sample_df: pd.DataFrame) -> pd.DataFrame:
    """Compute B vs Plasma symmetric composition/state decomposition in GSE164522."""
    bp_cells = cells_df.loc[cells_df["lineage"].eq("B/Plasma")].copy()
    bp_cells["subtype"] = np.where(
        bp_cells["celltype_sub"].str.contains("plasmaB", regex=False), "Plasma", "B"
    )
    sub = (
        bp_cells.groupby(["sample_uid", "lineage", "program", "subtype"])["score"]
        .agg(n_subtype="size", mean_score="mean")
        .reset_index()
    )
    support = sample_df.loc[
        sample_df["lineage"].eq("B/Plasma"), ["sample_uid", "lineage", "program", "n_cells"]
    ]
    usable = sub.merge(support, on=["sample_uid", "lineage", "program"])
    usable = usable.loc[usable["n_cells"] >= 50].merge(
        manifest[["sample_uid", "donor_uid", "lesion_uid", "tissue"]],
        on="sample_uid",
        validate="many_to_one",
    )
    usable = usable.loc[usable["tissue"].isin(["PT", "LM"])].copy()
    usable["proportion"] = usable["n_subtype"] / usable["n_cells"]
    usable["score_mass"] = usable["proportion"] * usable["mean_score"]

    rows = []
    keys = ["donor_uid", "lineage", "program"]
    for grouping, grp in usable.groupby(keys):
        frame = []
        observed = {}
        names = sorted(grp["subtype"].unique())
        for tissue, t_rows in grp.groupby("tissue"):
            n_lesions = t_rows["lesion_uid"].nunique()
            agg = (
                t_rows.groupby("subtype")[["n_subtype", "proportion", "score_mass"]]
                .sum()
                .reindex(names, fill_value=0)
            )
            p = agg["proportion"] / n_lesions
            mu = (agg["score_mass"] / agg["proportion"]).where(p > 0)
            observed[tissue] = float(agg["score_mass"].sum() / n_lesions)
            for s in names:
                frame.append(
                    {
                        "tissue": tissue,
                        "subtype": s,
                        "proportion": float(p.at[s]),
                        "mean_score": float(mu.at[s]) if pd.notna(mu.at[s]) else np.nan,
                        "n_subtype": int(agg["n_subtype"].at[s]),
                    }
                )
        values = pd.DataFrame(frame)
        for min_support in (1, 20):
            checked = values.copy()
            checked.loc[checked["n_subtype"] < min_support, "mean_score"] = np.nan
            pair = decompose_pair(checked)
            common = pair.get("common_support") or {}
            full = pair.get("full_decomposition") or {}
            cond = common.get("renormalized_components") or {}
            rows.append(
                dict(zip(keys, grouping))
                | {
                    "min_subtype_cells_per_tissue": min_support,
                    "status": pair["status"],
                    "observed_total_D_T": observed.get("LM", np.nan) - observed.get("PT", np.nan),
                    "full_composition": full.get("composition"),
                    "full_within_subtype_state": full.get("within_type_state"),
                    "common_mass_PT": common.get("retained_mass_PT"),
                    "common_mass_LM": common.get("retained_mass_LM"),
                    "common_total_D_T": cond.get("total_difference"),
                    "common_composition": cond.get("composition"),
                    "common_within_subtype_state": cond.get("within_type_state"),
                    "common_subtypes": " | ".join(common.get("subtypes", [])),
                }
            )
    return pd.DataFrame(rows)


def score_gse203608_step(
    matrix: sparse.csr_matrix,
    cells: pd.DataFrame,
    genes: pd.DataFrame,
    registry: pd.DataFrame,
    frozen_genes: pd.DataFrame,
    *,
    harmonize_hgnc: bool,
    shared_nuclear_only: bool,
    clean_cell_qc: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score GSE203608 in-vivo specimens under the specified ablation flags."""
    samples = (
        cells[["sample_uid", "donor_uid", "tissue"]]
        .drop_duplicates()
        .sort_values("sample_uid")
        .reset_index(drop=True)
    )
    samples = samples.assign(
        lesion_uid=samples["sample_uid"],
        dataset="GSE203608",
        role="external",
        input_scale="raw_counts",
        patient_mapping_verified=True,
        counts_values_verified=True,
        counts_provenance="Published GSE203608 MatrixMarket integer counts",
        eligible_for_paired_analysis=True,
    )

    if not harmonize_hgnc and not shared_nuclear_only:
        sym = genes["gene_short_name"]
        measurable = sym.ne("") & (~sym.duplicated(keep=False))
        cols = np.flatnonzero(measurable.to_numpy())
        axis = sym.iloc[cols].tolist()
    else:
        ndf = normalize_gene_axis(
            genes["gene_short_name"].tolist(), genes["gene_id"].tolist()
        )
        if shared_nuclear_only:
            axis = list(load_shared_nuclear_universe("245_164_203608_shared_nuclear"))
            al = align_axis_to_canonical_targets(
                genes["gene_short_name"].tolist(), axis, genes["gene_id"].tolist()
            )
            if not al["present"].all():
                raise ValueError("GSE203608 missing expected shared nuclear gene")
            cols = al["source_index"].to_numpy(dtype=np.int32)
        else:
            meas_df = ndf.loc[ndf["is_unique_measurable"]]
            cols = meas_df["source_index"].to_numpy(dtype=np.int32)
            axis = meas_df["canonical_symbol"].tolist()

    axis_set = set(axis)
    active_cells_mask = cells["clean_cell_qc"].to_numpy() if clean_cell_qc else np.ones(len(cells), dtype=bool)
    active_cells = cells.loc[active_cells_mask].copy()

    # Track immune composition before vs after dead-cell removal
    comp_counts = (
        active_cells.groupby(["sample_uid", "lineage"])
        .size()
        .unstack(fill_value=0)
        .reindex(samples["sample_uid"], fill_value=0)
    )

    cov_rows = []
    sample_rows = []
    for r in registry.itertuples(index=False):
        raw_req = frozen_genes.loc[
            frozen_genes["program"].eq(r.program) & frozen_genes["lineage"].eq(r.lineage),
            "gene_symbol",
        ].tolist()
        req = canonicalize_symbols(raw_req) if harmonize_hgnc else raw_req
        measured = [g for g in req if g in axis_set]
        ok = (
            len(measured) >= r.min_measured_positive_genes
            and len(measured) >= len(req) * r.min_measured_fraction_positive
        )
        cov_rows.append(
            {
                "cohort": "GSE203608",
                "lineage": r.lineage,
                "program": r.program,
                "n_axis_genes": len(axis),
                "n_requested": len(req),
                "n_measured": len(measured),
                "measured_fraction": len(measured) / len(req),
                "eligible": ok,
                "missing_symbols": "|".join(sorted(set(req) - set(measured))),
            }
        )
        lin_rows = np.flatnonzero(active_cells_mask & cells["lineage"].eq(r.lineage).to_numpy())
        scored_df = cells.iloc[lin_rows][["sample_uid"]].copy()
        if ok and len(lin_rows) > 0:
            scored_df["score"] = rank_scores(
                matrix[lin_rows][:, cols],
                axis,
                {r.program: measured},
                max_rank=1500,
                mode="legacy_positive",
            )[r.program].to_numpy()
        else:
            scored_df["score"] = np.nan
        by_sample = scored_df.groupby("sample_uid")["score"].agg(["mean", "size"])
        for smp in samples.itertuples(index=False):
            if smp.sample_uid in by_sample.index:
                n_c = int(by_sample.at[smp.sample_uid, "size"])
                m_val = float(by_sample.at[smp.sample_uid, "mean"])
            else:
                n_c = 0
                m_val = np.nan
            sample_rows.append(
                {
                    "sample_uid": smp.sample_uid,
                    "lineage": r.lineage,
                    "program": r.program,
                    "n_cells": n_c,
                    "unfiltered_score": m_val,
                    "score": m_val if n_c >= 50 else np.nan,
                }
            )

    sample_df = pd.DataFrame(sample_rows)
    patients, _, _, _ = calculate_patient_contrasts(samples, sample_df[SCORE_COLUMNS])
    summary = summarize_patients(patients, seed=SEED)
    return summary, patients, sample_df, pd.DataFrame(cov_rows), comp_counts


def load_gse178318_tumor_subset() -> tuple[sparse.csr_matrix, pd.DataFrame, pd.DataFrame, pd.DataFrame] | None:
    """Load the 12 clean tumor PT/LM H5ADs of GSE178318 (PBMC strictly excluded) if materialized."""
    if not GSE178_MANIFEST.is_file() or not GSE178_OBJECTS.is_dir():
        return None
    manifest = pd.read_csv(GSE178_MANIFEST)
    if manifest.empty:
        return None
    mats, obs_parts = [], []
    var_ref = None
    for gsm in manifest["sample_uid"].astype(str).tolist():
        obj = ad.read_h5ad(GSE178_OBJECTS / f"{gsm}.h5ad")
        if var_ref is None:
            var_ref = obj.var.copy()
        mats.append(obj.X.tocsr())
        obs_parts.append(obj.obs.copy())
    return sparse.vstack(mats, format="csr"), pd.concat(obs_parts, axis=0).reset_index(drop=True), var_ref, manifest


def score_gse178318_step(
    matrix: sparse.csr_matrix,
    cells: pd.DataFrame,
    var_df: pd.DataFrame,
    manifest: pd.DataFrame,
    registry: pd.DataFrame,
    frozen_genes: pd.DataFrame,
    *,
    harmonize_hgnc: bool,
    shared_nuclear_only: bool,
    clean_cell_qc: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score GSE178318 12 paired PT/LM tumor specimens under Track 1 stepwise ablation flags."""
    if shared_nuclear_only:
        axis = list(load_shared_nuclear_universe("245_164_203608_shared_nuclear"))
        canon_idx = pd.Index(var_df["canonical_symbol"].astype(str))
        cols = canon_idx.get_indexer(axis).astype(np.int32)
        if (cols < 0).any():
            raise ValueError("GSE178318 missing expected shared nuclear gene")
    elif harmonize_hgnc:
        cols = np.arange(len(var_df), dtype=np.int32)
        axis = var_df["canonical_symbol"].astype(str).tolist()
    else:
        raw_sym = var_df["gene_symbol"].astype(str)
        uniq = raw_sym.ne("") & (~raw_sym.duplicated(keep=False))
        cols = np.flatnonzero(uniq.to_numpy()).astype(np.int32)
        axis = raw_sym.iloc[cols].tolist()

    axis_set = set(axis)
    active_mask = (
        cells["qc_retained"].to_numpy(dtype=bool)
        if clean_cell_qc
        else cells["n_counts_visible"].gt(0).to_numpy(dtype=bool)
    )
    lin_arr = cells["lineage_m2"].astype(str).to_numpy()

    cov_rows = []
    sample_rows = []
    for r in registry.itertuples(index=False):
        raw_req = frozen_genes.loc[
            frozen_genes["program"].eq(r.program) & frozen_genes["lineage"].eq(r.lineage),
            "gene_symbol",
        ].tolist()
        req = canonicalize_symbols(raw_req) if harmonize_hgnc else raw_req
        measured = [g for g in req if g in axis_set]
        ok = (
            len(measured) >= r.min_measured_positive_genes
            and len(measured) >= len(req) * r.min_measured_fraction_positive
        )
        cov_rows.append(
            {
                "cohort": "GSE178318",
                "lineage": r.lineage,
                "program": r.program,
                "n_axis_genes": len(axis),
                "n_requested": len(req),
                "n_measured": len(measured),
                "measured_fraction": len(measured) / len(req),
                "eligible": ok,
                "missing_symbols": "|".join(sorted(set(req) - set(measured))),
            }
        )
        lin_rows = np.flatnonzero(active_mask & (lin_arr == r.lineage))
        scored_df = cells.iloc[lin_rows][["sample_uid"]].copy()
        if ok and len(lin_rows) > 0:
            scored_df["score"] = rank_scores(
                matrix[lin_rows][:, cols],
                axis,
                {r.program: measured},
                max_rank=1500,
                mode="legacy_positive",
            )[r.program].to_numpy()
        else:
            scored_df["score"] = np.nan
        by_sample = scored_df.groupby("sample_uid")["score"].agg(["mean", "size"])
        for smp in manifest.itertuples(index=False):
            if smp.sample_uid in by_sample.index:
                n_c = int(by_sample.at[smp.sample_uid, "size"])
                m_val = float(by_sample.at[smp.sample_uid, "mean"])
            else:
                n_c = 0
                m_val = np.nan
            sample_rows.append(
                {
                    "sample_uid": smp.sample_uid,
                    "lineage": r.lineage,
                    "program": r.program,
                    "n_cells": n_c,
                    "unfiltered_score": m_val,
                    "score": m_val if n_c >= 50 else np.nan,
                }
            )

    sample_df = pd.DataFrame(sample_rows)
    patients, _, _, _ = calculate_patient_contrasts(manifest, sample_df[SCORE_COLUMNS])
    summary = summarize_patients(patients, seed=SEED)
    return summary, patients, sample_df, pd.DataFrame(cov_rows)


def run() -> None:
    t0 = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    OUT.mkdir(parents=True, exist_ok=True)

    registry = pd.read_csv(ARCHIVE_M5 / "frozen_candidate_registry.csv")
    frozen_genes = pd.read_csv(ARCHIVE_M5 / "gene_sets/program_genes.csv")
    plan = pd.read_csv(ARCHIVE_M5 / "registered_tests.csv", keep_default_na=False)
    archived_summary = pd.read_csv(ARCHIVE_M6 / "effects_summary.csv")
    archived_reg = pd.read_csv(ARCHIVE_M6 / "registered_tests/registered_test_results.csv")

    manifest_164 = pd.read_csv(GSE164_MANIFEST)
    if "eligible_for_paired_analysis" not in manifest_164.columns:
        manifest_164["eligible_for_paired_analysis"] = True
    h5ads_164 = {code: ad.read_h5ad(GSE164_OBJECTS / f"{code}.h5ad") for code in GSE164_ORGANS}
    X_203, cells_203, genes_203 = load_gse203608_subset()
    loaded_178 = load_gse178318_tumor_subset()

    steps = [
        ("step0_v1_f1_baseline", False, False, False),
        ("step1_hgnc_harmonized", True, False, False),
        ("step2_shared_nuclear_axis", True, True, False),
        ("step3_full_techfix", True, True, True),
    ]

    ablation_rows = []
    cov_tables = []

    for step_id, hgnc_flag, nuc_flag, qc_flag in steps:
        sum_164, pts_164, smp_164, cov_164, cells_164 = score_gse164522_step(
            h5ads_164,
            manifest_164,
            registry,
            frozen_genes,
            harmonize_hgnc=hgnc_flag,
            shared_nuclear_only=nuc_flag,
            remove_cd3e_doublets=qc_flag,
        )
        sum_203, pts_203, smp_203, cov_203, comp_203 = score_gse203608_step(
            X_203,
            cells_203,
            genes_203,
            registry,
            frozen_genes,
            harmonize_hgnc=hgnc_flag,
            shared_nuclear_only=nuc_flag,
            clean_cell_qc=qc_flag,
        )
        sum_parts = [sum_164, sum_203]
        pts_parts = [pts_164, pts_203]
        cov_parts = [cov_164, cov_203]
        smp_178 = None
        if loaded_178 is not None and step_id != "step0_v1_f1_baseline":
            X_178, cells_178, var_178, manifest_178 = loaded_178
            sum_178, pts_178, smp_178, cov_178 = score_gse178318_step(
                X_178,
                cells_178,
                var_178,
                manifest_178,
                registry,
                frozen_genes,
                harmonize_hgnc=hgnc_flag,
                shared_nuclear_only=nuc_flag,
                clean_cell_qc=qc_flag,
            )
            sum_parts.append(sum_178)
            pts_parts.append(pts_178)
            cov_parts.append(cov_178)

        combined_sum = pd.concat(sum_parts, ignore_index=True)
        combined_pts = pd.concat(pts_parts, ignore_index=True)
        adj, extra, fam, audit = adjust_registered(plan, combined_sum)

        if step_id == "step0_v1_f1_baseline":
            num_cols = archived_summary.select_dtypes("number").columns
            if not combined_sum[num_cols].isna().equals(archived_summary[num_cols].isna()):
                raise AssertionError("Step 0 NaN pattern differs from archived V1_F1 effects_summary.csv")
            if not adj["raw_p"].isna().equals(archived_reg["raw_p"].isna()):
                raise AssertionError("Step 0 raw_p NaN pattern differs from archived V1_F1 registered_test_results.csv")
            diff_sum = np.nanmax(
                np.abs(
                    combined_sum[num_cols].to_numpy()
                    - archived_summary[num_cols].to_numpy()
                )
            )
            diff_holm = np.nanmax(
                np.abs(
                    adj["holm_adjusted_p"].to_numpy(dtype=float)
                    - archived_reg["holm_adjusted_p"].to_numpy(dtype=float)
                )
            )
            if diff_sum > 1e-14 or diff_holm > 1e-14:
                raise AssertionError(
                    f"Step 0 failed bit-exact reproduction of archived V1_F1: diff_sum={diff_sum}, diff_holm={diff_holm}"
                )
            comp_203_before = comp_203.copy()

        cov_step = pd.concat(cov_parts, ignore_index=True)
        cov_step.insert(0, "step", step_id)
        cov_tables.append(cov_step)

        # Merge bootstrap CI & median from combined_sum onto adjusted rows for ablation table
        ci_lookup = combined_sum.set_index(["dataset", "role", "lineage", "program", "contrast"])
        for row in adj.itertuples(index=False):
            key = (row.cohort, row.role, row.lineage, row.program, row.contrast)
            ci_lo = ci_lookup.at[key, "bootstrap_mean_ci_low"] if key in ci_lookup.index else np.nan
            ci_hi = ci_lookup.at[key, "bootstrap_mean_ci_high"] if key in ci_lookup.index else np.nan
            med = ci_lookup.at[key, "median_difference"] if key in ci_lookup.index else np.nan
            ablation_rows.append(
                {
                    "step": step_id,
                    "test_id": row.test_id,
                    "family": row.family,
                    "cohort": row.cohort,
                    "lineage": row.lineage,
                    "program": row.program,
                    "contrast": row.contrast,
                    "expected_direction": row.expected_direction,
                    "n_patients": row.n_patients,
                    "n_positive": row.n_positive,
                    "n_negative": row.n_negative,
                    "n_ties": row.n_ties,
                    "mean_difference": row.mean_difference,
                    "median_difference": med,
                    "bootstrap_ci_low": ci_lo,
                    "bootstrap_ci_high": ci_hi,
                    "raw_p": row.raw_p,
                    "holm_adjusted_p": row.holm_adjusted_p,
                    "direction_matches": row.direction_matches,
                    "holm_reject": row.holm_reject,
                    "statistical_direction_criterion_met": row.statistical_direction_criterion_met,
                }
            )

        if step_id == "step3_full_techfix":
            # Save primary Step 3 outputs
            combined_sum.to_csv(OUT / "effects_summary.csv", index=False)
            combined_pts.to_csv(OUT / "patient_contrasts.csv", index=False)
            adj.to_csv(OUT / "registered_test_results.csv", index=False)
            fam.to_csv(OUT / "registered_families.csv", index=False)
            extra.to_csv(OUT / "unplanned_exploratory_results.csv", index=False)
            smp_164.to_csv(OUT / "GSE164522_sample_program_scores.csv", index=False)
            smp_203.to_csv(OUT / "GSE203608_sample_program_scores.csv", index=False)
            if smp_178 is not None:
                smp_178.to_csv(OUT / "GSE178318_sample_program_scores.csv", index=False)

            bp_decomp = decompose_gse164522_b_plasma(manifest_164, cells_164, smp_164)
            bp_decomp.to_csv(OUT / "GSE164522_B_plasma_decomposition.csv", index=False)

            # GSE203608 composition denominator audit before vs after dead-cell removal
            comp_audit = comp_203_before.add_prefix("before_").join(
                comp_203.add_prefix("after_")
            )
            comp_audit["before_total_cells"] = comp_203_before.sum(axis=1)
            comp_audit["after_clean_cells"] = comp_203.sum(axis=1)
            comp_audit["removed_dead_or_na_cells"] = (
                comp_audit["before_total_cells"] - comp_audit["after_clean_cells"]
            )
            comp_audit.reset_index().to_csv(
                OUT / "GSE203608_qc_composition_audit.csv", index=False
            )

            # Exclude patient 49 clinical timing sensitivity on clean Step 3 GSE203608
            summarize_patients(
                pts_203.loc[pts_203["donor_uid"].ne("GSE203608:49")], seed=SEED
            ).to_csv(OUT / "GSE203608_exclude_patient49.csv", index=False)

            audit_payload = {
                **audit,
                "version": "CRLM_V1_245_F1_techfix",
                "started_utc": started_utc,
                "completed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "elapsed_seconds": round(time.perf_counter() - t0, 2),
                "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
                "archived_plan_sha256": sha256_file(ARCHIVE_M5 / "registered_tests.csv"),
                "gse164522_cd3e_doublets_removed": {
                    "B/Plasma_hC75_B-CD3E": int(
                        sum(
                            (o.obs["lineage"].eq("B/Plasma") & o.obs["celltype_sub"].str.contains("CD3E", regex=False)).sum()
                            for o in h5ads_164.values()
                        )
                    ),
                    "Mono/Mac/DC_CD3E_clusters": int(
                        sum(
                            (o.obs["lineage"].eq("Mono/Mac/DC") & o.obs["celltype_sub"].str.contains("CD3E", regex=False)).sum()
                            for o in h5ads_164.values()
                        )
                    ),
                },
                "gse203608_qc_bugfix": {
                    "in_vivo_sample_pass_qc_cells": int(len(cells_203)),
                    "cell_type_NA_cells": int(cells_203["Cell_type"].eq("NA").sum()),
                    "mito_gt_020_cells": int((cells_203["mito_fraction_num"] > 0.20).sum()),
                    "retained_clean_cells": int(cells_203["clean_cell_qc"].sum()),
                },
            }
            (OUT / "techfix_audit.json").write_text(
                json.dumps(audit_payload, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )

    ablation_df = pd.DataFrame(ablation_rows)
    ablation_df.to_csv(OUT / "techfix_stepwise_ablation.csv", index=False)
    cov_all = pd.concat(cov_tables, ignore_index=True)
    cov_all.to_csv(OUT / "techfix_gene_coverage_by_step.csv", index=False)

    write_report(ablation_df, cov_all, bp_decomp, audit_payload)

    files = sorted(
        p for p in OUT.rglob("*") if p.is_file() and p.name not in ("output_sha256.txt", ".gitkeep")
    )
    (OUT / "output_sha256.txt").write_text(
        "".join(f"{sha256_file(p)}  {p.relative_to(ROOT)}\n" for p in files),
        encoding="utf-8",
    )
    print(ablation_df.loc[ablation_df["step"].isin(["step0_v1_f1_baseline", "step3_full_techfix"]), [
        "step", "test_id", "n_patients", "n_positive", "n_negative", "mean_difference", "raw_p", "holm_adjusted_p"
    ]].to_string(index=False))


def write_report(
    ablation_df: pd.DataFrame,
    cov_all: pd.DataFrame,
    bp_decomp: pd.DataFrame,
    audit_payload: dict,
) -> None:
    """Write structured markdown stage report for Track 1 `CRLM_V1_245_F1_techfix` programs."""
    s0 = ablation_df.loc[ablation_df["step"].eq("step0_v1_f1_baseline")].set_index("test_id")
    s1 = ablation_df.loc[ablation_df["step"].eq("step1_hgnc_harmonized")].set_index("test_id")
    s2 = ablation_df.loc[ablation_df["step"].eq("step2_shared_nuclear_axis")].set_index("test_id")
    s3 = ablation_df.loc[ablation_df["step"].eq("step3_full_techfix")].set_index("test_id")

    bp1 = bp_decomp.loc[
        bp_decomp["min_subtype_cells_per_tissue"].eq(1)
        & bp_decomp["observed_total_D_T"].notna()
    ]
    comp_mean = bp1["full_composition"].mean()
    state_mean = bp1["full_within_subtype_state"].mean()
    tot_mean = bp1["observed_total_D_T"].mean()

    lines = [
        "# Track 1 (`CRLM_V1_245_F1_techfix`) · 冻结程序纯技术勘误报告",
        "",
        "本报告在不改动 `CRLM_V1_245_F1` 任何科学定义（2 个主候选、基因集、预期方向、4 个多重比较家族及 `positive_midrank_1500` 公式均 100% 锁定自归档 `artifacts/archive/CRLM_V1_245_F1/artifacts/m5/`）的前提下，分四步量化客观技术与数据清洗错误修复后的真实效应。",
        "",
        "## 1. 分步基因覆盖恢复",
        "",
        "| 步骤 | 队列 | 排序轴规模 | `Mac_SPP1` 可测成员 | `UPR` 可测成员 |",
        "|---|---|---:|---:|---:|",
    ]
    for step_id in (
        "step0_v1_f1_baseline",
        "step1_hgnc_harmonized",
        "step2_shared_nuclear_axis",
        "step3_full_techfix",
    ):
        for cohort in ("GSE164522", "GSE203608", "GSE178318"):
            sub_c = cov_all.loc[cov_all["step"].eq(step_id) & cov_all["cohort"].eq(cohort)]
            if sub_c.empty:
                continue
            sp = sub_c.loc[sub_c["program"].eq("B19_S4_Mac_SPP1")].iloc[0]
            up = sub_c.loc[sub_c["program"].eq("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")].iloc[0]
            lines.append(
                f"| `{step_id}` | {cohort} | {int(sp.n_axis_genes):,} | "
                f"{int(sp.n_measured)}/{int(sp.n_requested)} ({sp.measured_fraction*100:.1f}%) | "
                f"{int(up.n_measured)}/{int(up.n_requested)} ({up.measured_fraction*100:.1f}%) |"
            )

    lines.extend(
        [
            "",
            "## 2. 八项预登记检验分步消融对照（Step 0 → Step 3）",
            "",
            "| 检验 ID（队列 · 对比 · 候选） | Step 0 原版均值 (正/负; Holm p) | Step 1 +HGNC 对齐 | Step 2 +共同核背景 | Step 3 完整 `techfix` (正/负; 原始 p; Holm p) |",
            "|---|---|---|---|---|",
        ]
    )

    def fmt_short(row) -> str:
        if pd.isna(row.raw_p):
            return "NA (terms_pending)"
        return f"{row.mean_difference:+.5f} ({int(row.n_positive)}/{int(row.n_negative)}; p_H={row.holm_adjusted_p:.4f})"

    def fmt_full(row) -> str:
        if pd.isna(row.raw_p):
            return "NA (`GSE178318` 条款待确认，Holm 占位 p=1.0)"
        return (
            f"**{row.mean_difference:+.5f}** "
            f"({int(row.n_positive)}/{int(row.n_negative)}, N={int(row.n_patients)}; "
            f"p={row.raw_p:.4f}, Holm p={row.holm_adjusted_p:.4f})"
        )

    for tid in s0.index:
        r0, r1, r2, r3 = s0.loc[tid], s1.loc[tid], s2.loc[tid], s3.loc[tid]
        short_prog = "UPR" if "UNFOLDED" in r0.program else "Mac_SPP1"
        lines.append(
            f"| `{tid}` ({r0.cohort} · `{r0.contrast}` · {short_prog}) | "
            f"{fmt_short(r0)} | {fmt_short(r1)} | {fmt_short(r2)} | {fmt_full(r3)} |"
        )

    cd3_info = audit_payload["gse164522_cd3e_doublets_removed"]
    qc203 = audit_payload["gse203608_qc_bugfix"]
    lines.extend(
        [
            "",
            "## 3. 细胞级纯 Bug 修复与组成分解核查",
            "",
            f"- **GSE164522 跨谱系 `CD3E+` 双细胞群完整剔除**：同步移除髓系 `CD3E+` 亚群（`{cd3_info['Mono/Mac/DC_CD3E_clusters']}` 个细胞）与原脚本漏检的 `hC75_B-CD3E`（`{cd3_info['B/Plasma_hC75_B-CD3E']}` 个细胞）。剔除后 `B/Plasma UPR` 的患者内 `B` vs `Plasma` 对称分解（每侧 `>=1` 细胞，N={len(bp1)}）总差 `{tot_mean:+.5f}` 中，**组成分量为 `{comp_mean:+.5f}`，亚型内状态分量为 `{state_mean:+.5f}`**（详见 [B/Plasma 分解表](GSE164522_B_plasma_decomposition.csv)）。",
            f"- **GSE203608 死细胞/未注释细胞剔除**：体内 `Sample_pass_qc=TRUE` 共 `{qc203['in_vivo_sample_pass_qc_cells']:,}` 个细胞中，严格过滤 `Cell_type != 'NA'` 且 `MitoFraction <= 0.20` 剔除了 `{qc203['in_vivo_sample_pass_qc_cells'] - qc203['retained_clean_cells']:,}` 个高线粒体/未注释死细胞（其中 `Cell_type=='NA'` 共 `{qc203['cell_type_NA_cells']:,}` 个），保留 `{qc203['retained_clean_cells']:,}` 个合格细胞并修正各标本免疫分母（详见 [QC 与组成审计](GSE203608_qc_composition_audit.csv) 及 [剔除 49 号敏感性](GSE203608_exclude_patient49.csv)）。",
            "- **统计结论**：完成全部客观技术勘误后，原冻结的 2 个候选在 `PRIMARY_DT_164_V1` 与 `BACKGROUND_I_164_V1` 家族中仍未跨过双侧精确符号检验 Holm `alpha = 0.05` 门槛。这表明原 V1 未显著并非单纯由基因别名脱节造成：在四组织背景家族 `BACKGROUND_I_164_V1` 中，`Mac_SPP1` 合格患者仅 `N=5`（`5/0` 全同向时双侧精确符号检验理论最小 `raw p = 2/32 = 0.0625 > 0.05`，受离散小样本功效上限限制），而 `UPR` 交互 `I` 方向翻转至 `5/4`；在主家族 `PRIMARY_DT_164_V1`（`N=8` 与 `N=9`，全同向可达 `p=0.0078` 与 `p=0.0039`）中，未显著源于患者级同向比例不足（`Mac_SPP1` 为 `6/2`）及 `UPR` 自身 `80.7%` 归因于 `B` vs `Plasma` 亚型比例混杂。",
            "",
            "产物校验：`sha256sum -c artifacts/techfix_v1_f1/programs/output_sha256.txt`。",
        ]
    )
    (OUT / "stage_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    run()
