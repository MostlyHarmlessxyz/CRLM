"""Score frozen V1.1 definitions in the original in-vivo UPenn specimens (GSE203608).

Supports:
1. Official HGNC canonicalization + Ensembl-ID (`gene_id`) disambiguation of duplicate symbol
   pairs (`33,626 -> 33,658` unique measurable genes) and projection onto the 3-cohort shared
   nuclear universe `245_164_203608_shared_nuclear` (18,733 nuclear genes), preserving symbols
   already on the shared nuclear axis without context-free re-mapping.
2. Strict cell-level QC (`Cell_type != "NA"` and `MitoFraction <= 0.20` on the full raw cell
   annotation table), removing 12,136 dead / high-mitochondrial cells and recomputing dual
   immune denominators (`total_resolved_immune` and `core_three_lineage_immune`).
3. Frozen detection-matched background controls (`matched_controls.csv` inside `freeze_dir`)
   with strict `requires_debiased` enforcement (zero silent fallback or external refitting):
   - Primary metric: `debiased_ucell_scores` (`U'(P) - U'(P_ctrl)`).
   - Parallel sensitivity metrics: `log1p_cp10k_debiased` (counts normalized to CP10k on the
     18,733 shared nuclear axis), `ucell_raw`, and `legacy_positive` (`positive_midrank_1500`).
4. Clean author labels (`main`) PLUS unified two-stage marker-rule re-annotation (`annotation_m1.classify`),
   enabling direct evaluation of both broad lineages (`T/NK`, `Mono/Mac/DC`, `B/Plasma`) and fine
   subtypes (`CD4 T cell`, `CD8 T cell`, `NK cell`, `Macrophage`, `plasma cell`, `B cell`) at both
   `>=50` (primary) and `>=20` cell thresholds.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(ROOT / "starter") not in sys.path:
    sys.path.insert(1, str(ROOT / "starter"))

from crlm.annotation_m1 import classify  # noqa: E402
from crlm.hgnc_alias_map import (  # noqa: E402
    align_axis_to_canonical_targets,
    get_hgnc_mapper,
    load_shared_nuclear_universe,
    normalize_gene_axis,
)
from crlm.program_scores import (  # noqa: E402
    controls_from_dataframe,
    debiased_ucell_scores,
    legacy_positive_rank_scores,
    log1p_cp10k_scores,
    ucell_scores,
)
from paired_effects import calculate_patient_contrasts, summarize_patients  # noqa: E402

DATASET = "GSE203608"
RAW = ROOT / "data/raw" / DATASET
OUT = ROOT / "artifacts/m6" / DATASET
FROZEN = ROOT / "artifacts/m5"
AUDIT = ROOT / "starter/audit_data/late_2026_data"

SITES = {"Colon": "NC", "Liver": "NL", "LiverMet": "LM", "Tumor": "PT"}
MYELOID = {
    "Macrophage",
    "cDC2 (CD1C+)",
    "cDC1 (BATF3+)",
    "LAMP3+ DC (LAMP3+)",
    "pDC (LILRA4+)",
}
PARENT_LINEAGE_MAP = {
    "T/NK": "T/NK",
    "CD4 T cell": "T/NK",
    "CD8 T cell": "T/NK",
    "Double_neg_or_unsplit_T": "T/NK",
    "NK cell": "T/NK",
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
AUTHOR_RESOLVED_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast")
RULE_RESOLVED_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast")
CORE_THREE_IMMUNE = ("T/NK", "Mono/Mac/DC", "B/Plasma")
UNSPLIT_AUTHOR_T_SUBTYPES = ("CD4 T cell", "CD8 T cell", "NK cell", "Double_neg_or_unsplit_T")
FOCAL_SUBTYPES_BY_PARENT = {
    "T/NK": ("CD4 T cell", "CD8 T cell", "NK cell"),
    "Mono/Mac/DC": ("Macrophage", "Monocyte", "monocyte/macrophage"),
    "B/Plasma": ("plasma cell", "B cell"),
}
MAX_MITO_FRACTION = 0.20
SCORE_COLUMNS = ["sample_uid", "lineage", "program", "score", "n_cells"]
SEED = 20261121


def sha(path: Path) -> str:
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


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


def map_author_cells(cells: pd.DataFrame) -> pd.DataFrame:
    """Attach clean cell QC flags, broad lineage, and harmonized author fine subtype."""
    if "MitoFraction" not in cells.columns:
        raise KeyError("GSE203608 cell table missing required 'MitoFraction' column")
    out = cells.copy()
    mito_num = pd.to_numeric(out["MitoFraction"], errors="coerce")
    out["mito_fraction_num"] = mito_num
    out["pass_mito_qc"] = mito_num.notna() & (mito_num <= MAX_MITO_FRACTION)
    out["clean_cell_qc"] = out["Cell_type"].ne("NA") & out["pass_mito_qc"]

    broad = np.full(len(out), "Other/Unknown", dtype=object)
    ct = out["Cell_type"].astype(str)
    ct_mye = out["Cell_type_myeloid"].astype(str)
    is_epi = ct.str.startswith("Epithelial")
    is_nonimm = is_epi | ct.isin(["Fibroblast", "Endothelial", "Hepatocyte", "Smooth muscle", "Cholangiocyte"])
    broad[is_nonimm] = "Nonimmune"
    broad[ct.isin(["Plasma cell", "B cell"])] = "B/Plasma"
    broad[ct.eq("Myeloid") & ct_mye.isin(MYELOID)] = "Mono/Mac/DC"
    broad[ct.eq("T cell")] = "T/NK"
    broad[ct.eq("Mast cell")] = "Mast"
    out["lineage"] = pd.Series(broad, index=out.index, dtype=str)

    sub = np.full(len(out), "Other/Unknown", dtype=object)
    sub[is_epi] = "epithelial"
    sub[ct.eq("Fibroblast")] = "fibroblast"
    sub[ct.eq("Endothelial")] = "endothelial"
    sub[ct.eq("Hepatocyte")] = "hepatocyte"
    sub[ct.eq("Plasma cell")] = "plasma cell"
    sub[ct.eq("B cell")] = "B cell"
    sub[ct.eq("Myeloid") & ct_mye.eq("Macrophage")] = "Macrophage"
    sub[ct.eq("Myeloid") & ct_mye.isin(["cDC2 (CD1C+)", "cDC1 (BATF3+)"])] = "conventional DC"
    sub[ct.eq("Myeloid") & ct_mye.eq("LAMP3+ DC (LAMP3+)")] = "mature DC"
    sub[ct.eq("Myeloid") & ct_mye.eq("pDC (LILRA4+)")] = "plasmacytoid DC"
    sub[ct.eq("T cell")] = "T cell (author unsplit)"
    sub[ct.eq("Mast cell")] = "mast cell"
    out["author_subtype"] = pd.Series(sub, index=out.index, dtype=str)
    return out


def match_rule_population(
    target: str,
    rule_broad: np.ndarray,
    rule_sub: np.ndarray,
    rule_mf: np.ndarray,
) -> np.ndarray:
    if target in ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast", "Neutrophil", "Nonimmune"):
        return rule_broad == target
    if target == "monocyte/macrophage":
        return (rule_broad == "Mono/Mac/DC") & (rule_sub == "monocyte/macrophage")
    if target in ("Macrophage", "Monocyte"):
        return (rule_broad == "Mono/Mac/DC") & (rule_sub == "monocyte/macrophage") & (rule_mf == target)
    parent = PARENT_LINEAGE_MAP.get(target)
    if parent is not None:
        return (rule_broad == parent) & (rule_sub == target)
    return (rule_broad == target) | (rule_sub == target)


def match_author_population(
    target: str,
    auth_broad: np.ndarray,
    auth_sub: np.ndarray,
    rule_sub: np.ndarray | None = None,
    *,
    allow_stage2_t_split: bool = False,
) -> tuple[np.ndarray, str]:
    if target in ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast"):
        return auth_broad == target, "clean_author_exact"
    if target == "monocyte/macrophage":
        return (auth_broad == "Mono/Mac/DC") & (auth_sub == "Macrophage"), "clean_author_exact"
    if target in ("Macrophage", "conventional DC", "mature DC", "plasmacytoid DC", "B cell", "plasma cell", "mast cell"):
        return auth_sub == target, "clean_author_exact"
    if target in UNSPLIT_AUTHOR_T_SUBTYPES:
        if allow_stage2_t_split and rule_sub is not None:
            return (auth_broad == "T/NK") & (rule_sub == target), "author_T_cell_stage2_marker_split"
        return np.zeros(len(auth_broad), dtype=bool), "author_unsplit_T_unavailable"
    return np.zeros(len(auth_broad), dtype=bool), "unsupported_author_subtype"


def verify_freeze_integrity(freeze_dir: Path) -> int:
    """Verify SHA256 integrity of `freeze_dir` files (resolving relative paths strictly)."""
    freeze_resolved = Path(freeze_dir).resolve()
    root_resolved = ROOT.resolve()
    prog_scores_copy = freeze_resolved / "program_scores.py"
    if not prog_scores_copy.is_file():
        raise FileNotFoundError(f"Missing frozen program_scores.py in {freeze_resolved}")
    if sha(prog_scores_copy) != sha(root_resolved / "src/crlm/program_scores.py"):
        raise ValueError("Current program_scores.py differs from the frozen M5 copy")

    checked = 0
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
                in_freeze = True
            else:
                target = (root_resolved / rel_clean).resolve()
                in_freeze = target.is_relative_to(freeze_resolved)
            if in_freeze:
                if not target.is_file():
                    raise FileNotFoundError(f"Frozen M5 file declared in manifest is missing: {rel_clean} -> {target}")
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
        else set(load_shared_nuclear_universe("245_164_203608_shared_nuclear"))
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


def read_matrix(path: Path, nfeatures: int, ncells: int, gene_map: np.ndarray, cell_map: np.ndarray):
    blocks = []
    checked = 0
    with gzip.open(path, "rb") as f:
        if f.readline().strip() != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError("Expected author MatrixMarket integer count input")
        header = f.readline()
        while header.startswith(b"%"):
            header = f.readline()
        nf, nc, nnz = map(int, header.split())
        if (nf, nc) != (nfeatures, ncells):
            raise ValueError(f"Source annotation axes differ: {(nf, nc)} vs {(nfeatures, ncells)}")
        for chunk in pd.read_csv(
            f,
            sep=r"\s+",
            header=None,
            names=["g", "c", "count"],
            dtype=np.int64,
            chunksize=1_500_000,
        ):
            g = chunk.g.to_numpy() - 1
            c = chunk.c.to_numpy() - 1
            val = chunk["count"].to_numpy()
            if (
                g.min() < 0
                or g.max() >= nf
                or c.min() < 0
                or c.max() >= nc
                or val.min() < 0
                or val.max() > np.iinfo(np.int32).max
            ):
                raise ValueError("Invalid original count or coordinate")
            checked += len(val)
            use = (gene_map[g] >= 0) & (cell_map[c] >= 0) & (val > 0)
            if use.any():
                blocks.append((cell_map[c[use]], gene_map[g[use]], val[use].astype(np.int32)))
        if checked != nnz:
            raise ValueError(f"Read {checked} entries, source declares {nnz}")
    cell, gene, values = (np.concatenate(v) for v in zip(*blocks))
    matrix = sparse.coo_matrix(
        (values, (cell, gene)), shape=(int(cell_map.max()) + 1, int(gene_map.max()) + 1)
    ).tocsr()
    if matrix.nnz != len(values):
        raise ValueError("Duplicate original coordinates among selected cells and features")
    return matrix, nnz


def evaluate_in_vivo_cohort(
    X_full_canon: sparse.csr_matrix,
    full_canon_axis: Sequence[str],
    shared_axis: Sequence[str],
    in_vivo_cells: pd.DataFrame,
    samples: pd.DataFrame,
    registry: pd.DataFrame,
    sets: pd.DataFrame,
    controls_map: Mapping[tuple[str, str], Mapping[str, Sequence[str]]],
    *,
    unsplit_t_policy: str = "rule_within_author_T",
    requires_debiased: bool | None = None,
) -> dict[str, pd.DataFrame]:
    """Core scoring & contrast engine shared by production `run()` and synthetic `selfcheck`."""
    req_deb = registry_requires_debiased(registry) if requires_debiased is None else bool(requires_debiased)
    full_index = pd.Index([str(g) for g in full_canon_axis])
    shared_cols = full_index.get_indexer(shared_axis)
    if (shared_cols < 0).any():
        missing = [shared_axis[i] for i in np.flatnonzero(shared_cols < 0)[:10]]
        raise ValueError(f"GSE203608 canonical axis missing shared nuclear genes: {missing}")

    X_shared = X_full_canon[:, shared_cols].tocsr()
    total_counts = np.asarray(X_full_canon.sum(axis=1, dtype=np.int64)).ravel()
    nonempty = total_counts > 0

    rule_labels, rule_details = classify(X_full_canon, list(full_canon_axis))
    rule_broad = np.asarray(rule_labels, dtype=str)
    rule_sub = rule_details["subtype_hierarchical"].astype(str).to_numpy()
    rule_mf = rule_details["myeloid_fine"].fillna("").astype(str).to_numpy()

    cells = map_author_cells(in_vivo_cells).reset_index(drop=True)
    clean_author = cells["clean_cell_qc"].to_numpy(dtype=bool) & nonempty
    clean_mito_for_rule = clean_author
    auth_broad = cells["lineage"].astype(str).to_numpy()
    auth_sub = cells["author_subtype"].astype(str).to_numpy()
    sample_arr = cells["sample_uid"].astype(str).to_numpy()

    qc_rows = []
    comp_rows = []
    for samp in samples["sample_uid"].tolist():
        sm = sample_arr == samp
        n_in = int(sm.sum())
        n_na = int((sm & cells["Cell_type"].eq("NA").to_numpy()).sum())
        n_hi_mt = int((sm & (~cells["pass_mito_qc"].to_numpy(dtype=bool))).sum())
        n_cln = int((sm & clean_author).sum())
        n_rule_qc = int((sm & clean_mito_for_rule).sum())
        qc_rows.append(
            {
                "sample_uid": samp,
                "input_in_vivo_cells": n_in,
                "clean_author_cells": n_cln,
                "excluded_total_author": n_in - n_cln,
                "cells_author_NA": n_na,
                "cells_mito_gt_0p20_or_na": n_hi_mt,
                "clean_mito_cells_for_rule_m1": n_rule_qc,
            }
        )
        n_auth_tot = int((sm & clean_author & np.isin(auth_broad, AUTHOR_RESOLVED_IMMUNE)).sum())
        n_auth_c3 = int((sm & clean_author & np.isin(auth_broad, CORE_THREE_IMMUNE)).sum())
        n_rule_tot = int((sm & clean_mito_for_rule & np.isin(rule_broad, RULE_RESOLVED_IMMUNE)).sum())
        n_rule_c3 = int((sm & clean_mito_for_rule & np.isin(rule_broad, CORE_THREE_IMMUNE)).sum())
        for lin_name in ("T/NK", "Mono/Mac/DC", "B/Plasma", "Mast", "Neutrophil"):
            na_c = int((sm & clean_author & (auth_broad == lin_name)).sum())
            nr_c = int((sm & clean_mito_for_rule & (rule_broad == lin_name)).sum())
            comp_rows.append(
                {
                    "sample_uid": samp,
                    "lineage": lin_name,
                    "n_in_vivo_raw": n_in,
                    "n_clean_author_sample": n_cln,
                    "n_author_clean_cells": na_c,
                    "n_author_total_resolved_immune": n_auth_tot,
                    "n_author_core_three_immune": n_auth_c3,
                    "prop_author_total_resolved_immune": (
                        na_c / n_auth_tot if n_auth_tot > 0 and lin_name in AUTHOR_RESOLVED_IMMUNE else np.nan
                    ),
                    "prop_author_core_three_immune": (
                        na_c / n_auth_c3 if n_auth_c3 > 0 and lin_name in CORE_THREE_IMMUNE else np.nan
                    ),
                    "n_rule_m1_cells": nr_c,
                    "n_rule_total_resolved_immune": n_rule_tot,
                    "n_rule_core_three_immune": n_rule_c3,
                    "prop_rule_total_resolved_immune": (
                        nr_c / n_rule_tot if n_rule_tot > 0 and lin_name in RULE_RESOLVED_IMMUNE else np.nan
                    ),
                    "prop_rule_core_three_immune": (
                        nr_c / n_rule_c3 if n_rule_c3 > 0 and lin_name in CORE_THREE_IMMUNE else np.nan
                    ),
                }
            )

    concordance = (
        pd.DataFrame(
            {
                "clean_author_qc": clean_author,
                "pass_mito_qc": clean_mito_for_rule,
                "author_lineage": auth_broad,
                "author_subtype": auth_sub,
                "rule_lineage_m1": rule_broad,
                "rule_subtype_m1": rule_sub,
                "rule_myeloid_fine_m1": rule_mf,
            }
        )
        .value_counts()
        .rename("n_cells")
        .reset_index()
    )

    shared_set = set(shared_axis)
    measured_rows = []
    cell_score_frames = []

    for r in registry.itertuples(index=False):
        lin = str(r.lineage)
        prog = str(r.program)
        raw_req = sets.loc[sets.program.eq(prog) & sets.lineage.eq(lin), "gene_symbol"].tolist()
        canon_req = canonicalize_preserving_axis(raw_req, shared_set)
        keep = [g for g in canon_req if g in shared_set]
        ok = (
            len(keep) >= int(r.min_measured_positive_genes)
            and len(keep) >= len(canon_req) * float(r.min_measured_fraction_positive)
        )

        allow_split = unsplit_t_policy == "rule_within_author_T"
        m_main_raw, main_src = match_author_population(
            lin, auth_broad, auth_sub, rule_sub, allow_stage2_t_split=allow_split
        )
        m_strict_raw, _ = match_author_population(
            lin, auth_broad, auth_sub, rule_sub, allow_stage2_t_split=False
        )
        m_stage2_raw, _ = match_author_population(
            lin, auth_broad, auth_sub, rule_sub, allow_stage2_t_split=True
        )
        m_rule_raw = match_rule_population(lin, rule_broad, rule_sub, rule_mf)

        in_main = clean_author & m_main_raw
        in_strict = clean_author & m_strict_raw
        in_stage2 = clean_author & m_stage2_raw
        in_rule = clean_mito_for_rule & m_rule_raw
        in_uncleaned = nonempty & m_strict_raw

        measured_rows.append(
            dict(
                lineage=lin,
                program=prog,
                n_axis_genes=len(shared_axis),
                n_frozen_raw=len(raw_req),
                n_frozen=len(canon_req),
                n_measured=len(keep),
                measured_fraction=(len(keep) / len(canon_req) if canon_req else 0.0),
                has_frozen_matched_controls=(lin, prog) in controls_map,
                primary_annotation_source=main_src,
                n_cells_main=int(in_main.sum()),
                n_cells_rule_m1=int(in_rule.sum()),
                missing_symbols=" ".join(sorted(set(canon_req) - set(keep))),
                eligible=ok,
            )
        )

        if ok and req_deb and (lin, prog) not in controls_map:
            raise ValueError(
                f"Missing frozen matched_controls for {lin}::{prog}; silent fallback to non-debiased scores is prohibited"
            )

        union_rows = np.flatnonzero(in_main | in_strict | in_stage2 | in_rule | in_uncleaned)
        if not ok or len(union_rows) == 0:
            continue

        sub_X = X_shared[union_rows]
        leg_vec = legacy_positive_rank_scores(sub_X, shared_axis, {prog: keep}, max_rank=1500)[prog].to_numpy()
        if req_deb:
            c_spec = {prog: controls_map[(lin, prog)]}
            deb_df, raw_u_df, _ = debiased_ucell_scores(
                sub_X, shared_axis, {prog: keep}, matched_controls=c_spec, max_rank=1500, return_components=True
            )
            log_deb_vec = log1p_cp10k_scores(
                sub_X, shared_axis, {prog: keep}, matched_controls=c_spec, input_scale="counts"
            )[prog].to_numpy()
            deb_vec = deb_df[prog].to_numpy()
            raw_u_vec = raw_u_df[prog].to_numpy()
            primary_vec = deb_vec
        else:
            raw_u_vec = ucell_scores(sub_X, shared_axis, {prog: keep}, max_rank=1500)[prog].to_numpy()
            log_deb_vec = log1p_cp10k_scores(
                sub_X, shared_axis, {prog: keep}, debiased=False, input_scale="counts"
            )[prog].to_numpy()
            deb_vec = leg_vec
            primary_vec = leg_vec

        cell_score_frames.append(
            pd.DataFrame(
                {
                    "sample_uid": sample_arr[union_rows],
                    "donor_uid": cells["donor_uid"].to_numpy(str)[union_rows],
                    "tissue": cells["tissue"].to_numpy(str)[union_rows],
                    "lineage": lin,
                    "program": prog,
                    "author_subtype": auth_sub[union_rows],
                    "rule_subtype_m1": rule_sub[union_rows],
                    "rule_myeloid_fine_m1": rule_mf[union_rows],
                    "in_main": in_main[union_rows],
                    "in_author_strict": in_strict[union_rows],
                    "in_author_stage2_t": in_stage2[union_rows],
                    "in_rule_m1": in_rule[union_rows],
                    "in_uncleaned_author": in_uncleaned[union_rows],
                    "score": primary_vec,
                    "ucell_debiased": deb_vec,
                    "ucell_raw": raw_u_vec,
                    "log1p_cp10k_debiased": log_deb_vec,
                    "legacy_positive": leg_vec,
                }
            )
        )

    scored_cells = (
        pd.concat(cell_score_frames, ignore_index=True)
        if cell_score_frames
        else pd.DataFrame(
            columns=[
                "sample_uid", "donor_uid", "tissue", "lineage", "program",
                "author_subtype", "rule_subtype_m1", "rule_myeloid_fine_m1",
                "in_main", "in_author_strict", "in_author_stage2_t", "in_rule_m1",
                "in_uncleaned_author", "score", "ucell_debiased", "ucell_raw",
                "log1p_cp10k_debiased", "legacy_positive",
            ]
        )
    )

    def build_sample_table(mask_col: str, val_col: str) -> pd.DataFrame:
        sub = scored_cells.loc[scored_cells[mask_col]] if not scored_cells.empty else scored_cells
        agg = (
            sub.groupby(["sample_uid", "lineage", "program"])[val_col].agg(["mean", "size"]).reset_index()
            if not sub.empty
            else pd.DataFrame(columns=["sample_uid", "lineage", "program", "mean", "size"])
        )
        rows = []
        agg_idx = agg.set_index(["sample_uid", "lineage", "program"])
        for r in registry.itertuples(index=False):
            for s in samples.itertuples(index=False):
                key = (s.sample_uid, str(r.lineage), str(r.program))
                if key in agg_idx.index:
                    n_c = int(agg_idx.at[key, "size"])
                    val = float(agg_idx.at[key, "mean"]) if n_c > 0 else np.nan
                else:
                    n_c = 0
                    val = np.nan
                rows.append(
                    {
                        "sample_uid": s.sample_uid,
                        "lineage": str(r.lineage),
                        "program": str(r.program),
                        "n_cells": n_c,
                        "score": val,
                    }
                )
        return pd.DataFrame(rows)

    all_scores_main = build_sample_table("in_main", "score")

    policies = [
        ("main", "in_main", "score", 50),
        ("n20", "in_main", "score", 20),
        ("n100", "in_main", "score", 100),
        ("rule_annotation_m1", "in_rule_m1", "score", 50),
        ("rule_annotation_m1_n20", "in_rule_m1", "score", 20),
        ("rule_annotation_m1_log1p_debiased", "in_rule_m1", "log1p_cp10k_debiased", 50),
        ("log1p_cp10k_debiased", "in_main", "log1p_cp10k_debiased", 50),
        ("ucell_raw", "in_main", "ucell_raw", 50),
        ("legacy_positive", "in_main", "legacy_positive", 50),
        ("author_strict_only", "in_author_strict", "score", 50),
        ("author_T_stage2_split", "in_author_stage2_t", "score", 50),
        ("uncleaned_include_NA_highMT", "in_uncleaned_author", "score", 50),
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

    # Focal fine-subtype contrasts at both >=50 (primary) and >=20
    focal_summaries, focal_patients = [], []
    for min_cells in (50, 20):
        for branch_name, m_col, s_col, mf_col in (
            ("clean_author", "in_main", "author_subtype", None),
            ("rule_annotation_m1", "in_rule_m1", "rule_subtype_m1", "rule_myeloid_fine_m1"),
        ):
            b_sub = scored_cells.loc[scored_cells[m_col]] if not scored_cells.empty else scored_cells
            for r in registry.itertuples(index=False):
                lin, prog = str(r.lineage), str(r.program)
                parent_lin = PARENT_LINEAGE_MAP.get(lin, lin)
                targets = FOCAL_SUBTYPES_BY_PARENT.get(parent_lin, ())
                p_rows = b_sub.loc[b_sub["lineage"].eq(lin) & b_sub["program"].eq(prog)] if not b_sub.empty else b_sub
                for st in targets:
                    if p_rows.empty:
                        s_agg = pd.DataFrame(columns=["sample_uid", "mean", "size"])
                    else:
                        if st == "monocyte/macrophage":
                            smask = p_rows[s_col].isin(["monocyte/macrophage", "Macrophage", "Monocyte"])
                        elif st in ("Macrophage", "Monocyte") and mf_col is not None:
                            smask = p_rows[s_col].eq(st) | p_rows[mf_col].eq(st)
                        else:
                            smask = p_rows[s_col].eq(st)
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
                    p_df.insert(0, "annotation_branch", branch_name)
                    s_df.insert(0, "min_cells_threshold", min_cells)
                    s_df.insert(0, "parent_candidate_lineage", lin)
                    s_df.insert(0, "annotation_branch", branch_name)
                    focal_patients.append(p_df)
                    focal_summaries.append(s_df)

    captured_clean = (
        cells.loc[clean_author]
        .groupby(["sample_uid", "lineage"])
        .size()
        .unstack(fill_value=0)
        .reindex(samples["sample_uid"], fill_value=0)
    )

    return {
        "cell_qc_audit": pd.DataFrame(qc_rows),
        "sample_immune_composition": pd.DataFrame(comp_rows),
        "captured_composition_counts": captured_clean,
        "concordance": concordance,
        "measurement_coverage": pd.DataFrame(measured_rows),
        "all_scores_main": all_scores_main,
        "sensitivity_summary": pd.concat(sens_summaries, ignore_index=True),
        "focal_subtype_effects_summary": (
            pd.concat(focal_summaries, ignore_index=True) if focal_summaries else pd.DataFrame()
        ),
        "focal_subtype_patient_contrasts": (
            pd.concat(focal_patients, ignore_index=True) if focal_patients else pd.DataFrame()
        ),
        "policy_outputs": policy_outputs,
        "n_clean_author_cells": int(clean_author.sum()),
        "n_excluded_author_cells": int(len(cells) - int(clean_author.sum())),
    }


def load_verified_cell_table(cell_file: Path, subset_file: Path, expected_sha256: str) -> pd.DataFrame:
    """Verify full raw cell annotation table against archived identity subset while preserving MitoFraction."""
    if sha(cell_file) != expected_sha256:
        raise ValueError("Original author annotation version differs from identity audit")
    raw_cells = pd.read_csv(cell_file, dtype=str, keep_default_na=False)
    identity_subset = pd.read_csv(subset_file, dtype=str, keep_default_na=False)
    cells = raw_cells.rename(columns={raw_cells.columns[0]: "barcode"})
    pd.testing.assert_frame_equal(cells[identity_subset.columns], identity_subset)
    if not cells.barcode.is_unique:
        raise ValueError("Duplicate author barcode")
    if "MitoFraction" not in cells.columns:
        raise KeyError("Full GSE203608 cell annotation table missing 'MitoFraction'")
    return cells


def run(
    freeze_dir: Path = FROZEN,
    out_dir: Path = OUT,
    unsplit_t_policy: str = "rule_within_author_T",
):
    start = time.perf_counter()
    freeze_dir = Path(freeze_dir).resolve()
    out_dir = Path(out_dir).resolve()

    verify_freeze_integrity(freeze_dir)
    registry = pd.read_csv(freeze_dir / "frozen_candidate_registry.csv")
    if registry.empty or not registry.status.astype(str).str.lower().eq("frozen").all():
        raise ValueError("A real frozen registry is required")
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
            extra_pairs.extend(zip(sec_df.lineage.astype(str), sec_df.program.astype(str)))
        extra_pairs.append(("B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"))
        mg_keys = set(zip(mg_df.lineage.astype(str), mg_df.program.astype(str)))
        existing = set(zip(eval_registry.lineage.astype(str), eval_registry.program.astype(str)))
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

    shared_axis = list(load_shared_nuclear_universe("245_164_203608_shared_nuclear"))
    cand_keys = sorted(set(zip(eval_registry.lineage.astype(str), eval_registry.program.astype(str))))
    controls_map = load_frozen_matched_controls(
        freeze_dir, cand_keys, registry, shared_axis=shared_axis, strict_keys=cand_keys
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    coverage_dir = out_dir / "coverage"
    coverage_dir.mkdir(exist_ok=True)

    cell_file = RAW / "GSE203608_cell_annotation.csv.gz"
    subset_file = AUDIT / "GSE203608_original_identity_annotation_subset.csv.gz"
    source_meta = json.loads((AUDIT / "GSE203608_download_provenance.json").read_text())[0]
    cells = load_verified_cell_table(cell_file, subset_file, source_meta["sha256"])

    in_vivo = cells.SampleType.isin(SITES) & cells.Sample_pass_qc.eq("TRUE")
    in_vivo_cells = cells.loc[in_vivo].copy()
    in_vivo_cells["sample_uid"] = DATASET + ":" + in_vivo_cells.Dataset
    in_vivo_cells["donor_uid"] = DATASET + ":" + in_vivo_cells.Patient
    in_vivo_cells["tissue"] = in_vivo_cells.SampleType.map(SITES)

    samples = in_vivo_cells[["sample_uid", "donor_uid", "tissue"]].drop_duplicates()
    if samples.sample_uid.duplicated().any():
        raise ValueError("A specimen maps to more than one patient or tissue")
    samples = samples.assign(
        lesion_uid=samples.sample_uid,
        dataset=DATASET,
        role="external",
        input_scale="raw_counts",
        patient_mapping_verified=True,
        counts_values_verified=True,
        counts_provenance="Published original gene_by_cell_count_matrix and ordered annotations; checked integer entries",
    )

    gene_file = RAW / "GSE203608_gene_annotation.csv.gz"
    genes = pd.read_csv(gene_file, dtype=str, keep_default_na=False)
    if genes.gene_id.duplicated().any():
        raise ValueError("Duplicate author gene IDs")

    norm_genes = normalize_gene_axis(
        genes["gene_short_name"].tolist(),
        genes["gene_id"].tolist(),
        excluded_canonical=set(load_shared_nuclear_universe("excluded_353_canonical")),
    )
    measurable_mask = norm_genes["is_unique_measurable"].to_numpy(dtype=bool)
    gene_map = np.full(len(genes), -1, dtype=np.int32)
    gene_map[measurable_mask] = np.arange(int(measurable_mask.sum()), dtype=np.int32)
    full_canon_axis = norm_genes.loc[measurable_mask, "canonical_symbol"].tolist()

    al_shared = align_axis_to_canonical_targets(
        genes["gene_short_name"].tolist(), shared_axis, genes["gene_id"].tolist()
    )
    if not al_shared["present"].all():
        raise ValueError("GSE203608 missing expected genes from 245_164_203608_shared_nuclear")

    in_vivo_arr = in_vivo.to_numpy(dtype=bool)
    cell_map = np.full(len(cells), -1, dtype=np.int32)
    cell_map[in_vivo_arr] = np.arange(int(in_vivo_arr.sum()), dtype=np.int32)
    matrix_file = RAW / "GSE203608_gene_by_cell_count_matrix.txt.gz"
    X_full, nnz = read_matrix(matrix_file, len(genes), len(cells), gene_map, cell_map)

    res = evaluate_in_vivo_cohort(
        X_full,
        full_canon_axis,
        shared_axis,
        in_vivo_cells,
        samples,
        eval_registry,
        sets,
        controls_map,
        unsplit_t_policy=unsplit_t_policy,
        requires_debiased=registry_requires_debiased(registry),
    )

    samples.to_csv(coverage_dir / "sample_manifest.csv", index=False)
    res["captured_composition_counts"].to_csv(coverage_dir / "captured_composition_counts.csv")
    res["cell_qc_audit"].to_csv(coverage_dir / "cell_qc_audit.csv", index=False)
    res["sample_immune_composition"].to_csv(coverage_dir / "sample_immune_composition.csv", index=False)
    res["concordance"].to_csv(coverage_dir / "author_vs_rule_annotation_concordance.csv", index=False)
    res["measurement_coverage"].to_csv(coverage_dir / "measurement_coverage.csv", index=False)
    (coverage_dir / "external_rank_background.txt").write_text("\n".join(shared_axis) + "\n")
    res["all_scores_main"].to_csv(out_dir / "sample_program_measurements.csv", index=False)

    folder_map = {
        "main": "patient_effects",
        "n20": "sensitivity_20",
        "n100": "sensitivity_100",
        "rule_annotation_m1": "sensitivity_rule_annotation_m1",
        "log1p_cp10k_debiased": "sensitivity_log1p_cp10k_debiased",
        "legacy_positive": "sensitivity_legacy_positive",
    }
    for pol_name, folder_name in folder_map.items():
        p_out = res["policy_outputs"][pol_name]
        sub_dir = out_dir / folder_name
        sub_dir.mkdir(exist_ok=True)
        for fname in ("patient_contrasts", "effects_summary", "patient_tissue_means", "lesion_means"):
            p_out[fname].to_csv(sub_dir / f"{fname}.csv", index=False)
        if pol_name == "main":
            pat_main = p_out["patient_contrasts"]
            summarize_patients(pat_main[pat_main.donor_uid.ne(DATASET + ":49")], seed=SEED).to_csv(
                sub_dir / "exclude_patient49.csv", index=False
            )
            res["sensitivity_summary"].to_csv(sub_dir / "sensitivity_summary.csv", index=False)
            if not res["focal_subtype_effects_summary"].empty:
                res["focal_subtype_effects_summary"].to_csv(
                    sub_dir / "focal_subtype_effects_summary.csv", index=False
                )
                res["focal_subtype_patient_contrasts"].to_csv(
                    sub_dir / "focal_subtype_patient_contrasts.csv", index=False
                )

    profile = dict(
        completed_utc=utc_now(),
        source_metadata_file=cell_file.name,
        metadata_sha256=source_meta["sha256"],
        gene_annotation_sha256=sha(gene_file),
        expression_sha256=sha(matrix_file),
        code_sha256=sha(Path(__file__)),
        rank_code_sha256=sha(ROOT / "src/crlm/program_scores.py"),
        hgnc_code_sha256=sha(ROOT / "src/crlm/hgnc_alias_map.py"),
        annotation_code_sha256=sha(ROOT / "src/crlm/annotation_m1.py"),
        source_coordinates=nnz,
        source_rows=len(genes),
        source_cells=len(cells),
        unique_measurable_canonical_genes=len(full_canon_axis),
        shared_nuclear_axis_genes=len(shared_axis),
        in_vivo_cells_raw=int(in_vivo.sum()),
        in_vivo_cells_clean_author=res["n_clean_author_cells"],
        in_vivo_cells_excluded_NA_or_highMT=res["n_excluded_author_cells"],
        in_vivo_specimens=len(samples),
        donors=int(samples.donor_uid.nunique()),
        elapsed_seconds=round(time.perf_counter() - start, 3),
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
    )
    (out_dir / "resource_profile.json").write_text(json.dumps(profile, indent=2) + "\n")
    print(json.dumps(profile), flush=True)


def run_selfcheck() -> None:
    """Pure in-memory synthetic & real-metadata schema verification (zero project disk writes)."""
    import os
    import tempfile
    from crlm.annotation_m1 import AUXILIARY_MARKERS, PANELS

    rng = np.random.default_rng(SEED)
    full_shared = list(load_shared_nuclear_universe("245_164_203608_shared_nuclear"))
    full_shared_set = set(full_shared)

    # 1. Verify axis-preserving canonicalization keeps PRR34-AS1 / LINC01358 intact
    for tricky in ("PRR34-AS1", "LINC01358", "APOE", "APOC1"):
        if tricky in full_shared_set:
            assert canonicalize_preserving_axis([tricky], full_shared_set) == [tricky]

    # 2. Verify real `load_verified_cell_table` contract (preserving `MitoFraction` & `Epithelial(...)` mapping)
    cell_file = RAW / "GSE203608_cell_annotation.csv.gz"
    subset_file = AUDIT / "GSE203608_original_identity_annotation_subset.csv.gz"
    if cell_file.is_file() and subset_file.is_file():
        source_meta = json.loads((AUDIT / "GSE203608_download_provenance.json").read_text())[0]
        verified_cells = load_verified_cell_table(cell_file, subset_file, source_meta["sha256"])
        assert "MitoFraction" in verified_cells.columns and len(verified_cells) == 126014
        mapped_real = map_author_cells(verified_cells.loc[verified_cells.SampleType.isin(SITES) & verified_cells.Sample_pass_qc.eq("TRUE")])
        assert int(mapped_real["clean_cell_qc"].sum()) == 47525
        epi_rows = mapped_real.loc[mapped_real["Cell_type"].str.startswith("Epithelial")]
        assert not epi_rows.empty and epi_rows["lineage"].eq("Nonimmune").all()

    # 3. Verify freeze_dir resolution & missing-control ValueError guard
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
        reg_two = pd.DataFrame(
            [
                {"lineage": "Mono/Mac/DC", "program": "SYN_SAM", "scoring_method": "ucell_debiased"},
                {"lineage": "T/NK", "program": "MISSING_PROG", "scoring_method": "ucell_debiased"},
            ]
        )
        raised_val = False
        try:
            load_frozen_matched_controls(
                tdir, [("Mono/Mac/DC", "SYN_SAM"), ("T/NK", "MISSING_PROG")], reg_two, shared_axis=full_shared
            )
        except ValueError:
            raised_val = True
        assert raised_val, "Expected ValueError when a candidate program is missing from matched_controls.csv"

    shared_axis = full_shared[:220]
    marker_inject = [
        "CD3D", "CD3E", "CD4", "CD8A", "TRAC", "IL7R",
        "MS4A1", "CD79A", "MZB1", "JCHAIN", "SDC1", "TNFRSF17",
        "C1QA", "C1QB", "CD14", "CD68", "CSF1R", "SPP1", "PTPRC", "EPCAM", "KRT19",
    ]
    for i, mg in enumerate(marker_inject):
        if mg not in shared_axis:
            shared_axis[i] = mg
    all_panel_genes = {g for _, a, s in PANELS.values() for g in (a + " " + s).split()} | set(AUXILIARY_MARKERS)
    non_panel = [g for g in shared_axis if g not in all_panel_genes]
    prog_genes = non_panel[:10]
    ctrl_pool = non_panel[10:120]
    controls_map = {
        ("Mono/Mac/DC", "SYN_SAM"): {
            g: [ctrl_pool[(idx * 10 + r) % len(ctrl_pool)] for r in range(10)]
            for idx, g in enumerate(prog_genes)
        },
        ("CD8 T cell", "SYN_CD8"): {
            g: [ctrl_pool[(idx * 10 + r + 3) % len(ctrl_pool)] for r in range(10)]
            for idx, g in enumerate(prog_genes)
        },
    }
    registry = pd.DataFrame(
        [
            {
                "lineage": "Mono/Mac/DC",
                "program": "SYN_SAM",
                "status": "frozen",
                "scoring_method": "ucell_debiased",
                "min_measured_positive_genes": 5,
                "min_measured_fraction_positive": 0.8,
            },
            {
                "lineage": "CD8 T cell",
                "program": "SYN_CD8",
                "status": "frozen",
                "scoring_method": "ucell_debiased",
                "min_measured_positive_genes": 5,
                "min_measured_fraction_positive": 0.8,
            },
        ]
    )
    sets = pd.DataFrame(
        [{"lineage": "Mono/Mac/DC", "program": "SYN_SAM", "gene_symbol": g} for g in prog_genes]
        + [{"lineage": "CD8 T cell", "program": "SYN_CD8", "gene_symbol": g} for g in prog_genes]
    )

    col_idx = {g: i for i, g in enumerate(shared_axis)}
    mats, cell_frames, samp_rows = [], [], []
    for donor in ("10", "20"):
        for site, tissue in (("Tumor", "PT"), ("LiverMet", "LM")):
            s_uid = f"{DATASET}:P{donor}_{tissue}"
            samp_rows.append(
                {
                    "sample_uid": s_uid,
                    "donor_uid": f"{DATASET}:{donor}",
                    "lesion_uid": s_uid,
                    "tissue": tissue,
                    "dataset": DATASET,
                    "role": "external",
                    "input_scale": "raw_counts",
                    "patient_mapping_verified": True,
                }
            )
            n_c = 160  # 60 Macrophage + 60 unsplit T cell (CD8A+) + 10 Epithelial(Tumor) + 30 dead NA/high-MT
            m = sparse.random(n_c, len(shared_axis), density=0.15, format="csr", random_state=rng)
            m.data = rng.integers(1, 12, size=m.nnz).astype(np.int32)
            ml = m.tolil()
            for g in all_panel_genes:
                if g in col_idx:
                    ml[:, col_idx[g]] = 0
            for r in range(60):
                for g in ("C1QA", "C1QB", "CD68", "CSF1R"):
                    ml[r, col_idx[g]] = 15
                if tissue == "LM":
                    for g in prog_genes[:6]:
                        ml[r, col_idx[g]] = 30
            for r in range(60, 120):
                for g in ("CD3D", "CD3E", "CD8A", "TRAC"):
                    ml[r, col_idx[g]] = 15
                if tissue == "LM":
                    for g in prog_genes[:6]:
                        ml[r, col_idx[g]] = 25
            for r in range(120, 130):
                for g in ("EPCAM", "KRT19"):
                    ml[r, col_idx[g]] = 15
            mats.append(ml.tocsr())
            cell_frames.append(
                pd.DataFrame(
                    {
                        "sample_uid": s_uid,
                        "donor_uid": f"{DATASET}:{donor}",
                        "tissue": tissue,
                        "SampleType": site,
                        "Sample_pass_qc": "TRUE",
                        "Cell_type": ["Myeloid"] * 60 + ["T cell"] * 60 + ["Epithelial(Tumor)"] * 10 + ["NA"] * 30,
                        "Cell_type_myeloid": ["Macrophage"] * 60 + ["NA"] * 100,
                        "MitoFraction": ["0.05"] * 130 + ["0.45"] * 30,
                    }
                )
            )

    X_all = sparse.vstack(mats, format="csr")
    cells_all = pd.concat(cell_frames, ignore_index=True)
    samples_df = pd.DataFrame(samp_rows)
    res = evaluate_in_vivo_cohort(
        X_all, shared_axis, shared_axis, cells_all, samples_df, registry, sets, controls_map, requires_debiased=True
    )
    assert res["n_excluded_author_cells"] == 120 and res["n_clean_author_cells"] == 520
    # Verify both >=50 and >=20 focal subtype rows exist
    assert set(res["focal_subtype_effects_summary"]["min_cells_threshold"].unique()) == {20, 50}
    sens = res["sensitivity_summary"]
    cd8_rule = sens.loc[
        sens["analysis"].eq("rule_annotation_m1") & sens["lineage"].eq("CD8 T cell") & sens["contrast"].eq("D_T")
    ].iloc[0]
    cd8_strict = sens.loc[
        sens["analysis"].eq("author_strict_only") & sens["lineage"].eq("CD8 T cell") & sens["contrast"].eq("D_T")
    ].iloc[0]
    assert int(cd8_rule["n_patients"]) == 2 and float(cd8_rule["mean_difference"]) > 0
    assert int(cd8_strict["n_patients"]) == 0 and pd.isna(cd8_strict["sign_p_two_sided"])
    for _, grp in res["sample_immune_composition"].groupby("sample_uid"):
        assert np.isclose(grp["prop_author_total_resolved_immune"].dropna().sum(), 1.0)
        assert np.isclose(grp["prop_author_core_three_immune"].dropna().sum(), 1.0)
        assert np.isclose(grp["prop_rule_total_resolved_immune"].dropna().sum(), 1.0)
        assert np.isclose(grp["prop_rule_core_three_immune"].dropna().sum(), 1.0)
    print(
        json.dumps(
            {
                "selfcheck": "passed",
                "script": "validate_gse203608.py",
                "real_cell_table_verified": cell_file.is_file(),
                "clean_cells_synthetic": res["n_clean_author_cells"],
                "excluded_dead_cells_synthetic": res["n_excluded_author_cells"],
                "cd8_rule_m1_mean_DT": round(float(cd8_rule["mean_difference"]), 6),
            }
        ),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selfcheck", action="store_true", help="Run pure in-memory synthetic & schema contract check.")
    parser.add_argument("--freeze-dir", type=Path, default=FROZEN)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    parser.add_argument(
        "--unsplit-t-policy",
        choices=("rule_within_author_T", "strict_na"),
        default="rule_within_author_T",
    )
    args = parser.parse_args()
    if args.selfcheck:
        run_selfcheck()
    else:
        run(
            freeze_dir=args.freeze_dir.resolve(),
            out_dir=args.out_dir.resolve(),
            unsplit_t_policy=args.unsplit_t_policy,
        )


if __name__ == "__main__":
    main()
