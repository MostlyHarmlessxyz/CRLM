"""Deterministic tumor specimen bag builders and full-specimen composition helpers.

Guarantees:
1. Only verified tumor libraries (`PT` or `LM`) can form classification bags;
   normal tissues (`NC`, `NL`) raise `ValueError`.
2. By default (`max_cells=None`), retains ALL `qc_retained` cells for each eligible
   population (`population >= min_cells`), eliminating single-draw 256-cell truncation.
3. Supports deterministic `max_cells` subsampling with `repeat_index` and in-memory
   `subsample_bag_indices` for multi-draw (`256 x 10`) subsampling variance audits.
4. Supports `All3Lineages` joint multi-lineage specimen bags with clean hierarchical
   fine cell-type labels (`fine_cell_type`) for specimen-level HA/CTA/MixMIL.
5. Slices onto canonical HGNC symbols or a specified shared nuclear gene universe
   (e.g. `245_164_shared_nuclear`, 18,815 genes) with strict uniqueness and presence checks.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd

from crlm.hgnc_alias_map import (
    align_axis_to_canonical_targets,
    load_shared_nuclear_universe,
)

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261121
LINEAGES = ("T/NK", "Mono/Mac/DC", "B/Plasma")
JOINT_LINEAGE = "All3Lineages"
RESOLVED_IMMUNE_LINEAGES = ("T/NK", "Mono/Mac/DC", "B/Plasma", "Neutrophil", "Mast")
TUMOR_TISSUES = ("PT", "LM")
INVALID_SUBTYPE_LABELS = frozenset({"", "unresolved", "Unknown", "Cross_lineage_doublet", "nan", "None"})


def cell_seed(uid: object, seed: int = SEED) -> int:
    """Stable 32-bit random seed for a cell or bag, independent of processing order."""
    return int.from_bytes(
        hashlib.sha256(f"{seed}|{uid}".encode("utf-8")).digest()[:4], "little"
    )


@dataclass
class Bag:
    obs: pd.DataFrame
    counts: object
    genes: pd.Index
    population: int


def select_bag_rows(
    cell_ids: pd.Index,
    candidate_rows: np.ndarray,
    bag_key: str,
    *,
    min_cells: int = 50,
    max_cells: int | None = None,
    repeat_index: int = 0,
    seed: int = SEED,
) -> np.ndarray:
    """Deterministically select row indices from `candidate_rows`."""
    if min_cells < 1:
        raise ValueError("min_cells must be >= 1")
    if max_cells is not None and max_cells < 1:
        raise ValueError("max_cells must be None or >= 1")
    rows = np.asarray(candidate_rows, dtype=np.int64)
    if len(rows) < min_cells:
        return np.array([], dtype=np.int64)
    if max_cells is None or len(rows) <= max_cells:
        return np.sort(rows)
    ordered = rows[np.argsort(cell_ids[rows].to_numpy(dtype=str), kind="mergesort")]
    seed_key = bag_key if int(repeat_index) == 0 else f"{bag_key}|rep{int(repeat_index)}"
    rng = np.random.default_rng(cell_seed(seed_key, seed=seed))
    return np.sort(rng.choice(ordered, int(max_cells), replace=False))


def subsample_bag_indices(
    cell_uids: Sequence[str] | pd.Index,
    sample_uid: str,
    lineage: str,
    *,
    max_cells: int = 256,
    n_repeats: int = 10,
    seed: int = SEED,
) -> list[np.ndarray]:
    """Return `n_repeats` deterministic sub-bag index arrays (0-based into `cell_uids`).

    Repeat `r=0` matches `select_bag_rows(..., max_cells=max_cells, repeat_index=0)`
    on the same `cell_uids` universe.
    """
    if max_cells < 1 or n_repeats < 1:
        raise ValueError("max_cells and n_repeats must be >= 1")
    uid_index = pd.Index([str(u) for u in cell_uids])
    if len(uid_index) == 0 or not uid_index.is_unique:
        raise ValueError("cell_uids must be nonempty and unique")
    all_pos = np.arange(len(uid_index), dtype=np.int64)
    key = f"{sample_uid}|{lineage}"
    return [
        select_bag_rows(
            uid_index,
            all_pos,
            key,
            min_cells=1,
            max_cells=max_cells,
            repeat_index=r,
            seed=seed,
        )
        for r in range(n_repeats)
    ]


def _resolve_gene_slice(
    var: pd.DataFrame,
    *,
    use_canonical_genes: bool = True,
    gene_universe: str | Sequence[str] | None = None,
    source_label: str = "",
) -> tuple[np.ndarray | None, pd.Index]:
    has_canon = use_canonical_genes and ("canonical_symbol" in var.columns)
    base_col = "canonical_symbol" if has_canon else "gene_symbol"
    base_syms = pd.Index(var[base_col].astype(str))
    if not base_syms.is_unique or (base_syms == "").any():
        raise ValueError(f"Non-unique or empty {base_col} axis in {source_label}")

    if gene_universe is None:
        return None, base_syms

    targets = (
        list(load_shared_nuclear_universe(gene_universe))
        if isinstance(gene_universe, str)
        else [str(g) for g in gene_universe]
    )
    target_idx = pd.Index(targets)
    if len(target_idx) == 0 or not target_idx.is_unique:
        raise ValueError("Requested gene_universe must be nonempty and unique")

    if has_canon:
        col_indices = base_syms.get_indexer(target_idx)
        if (col_indices < 0).any():
            missing = target_idx[col_indices < 0][:10].tolist()
            raise ValueError(f"Missing {len(target_idx[col_indices < 0])} target genes in {source_label}: {missing}")
        return col_indices.astype(np.int64), target_idx

    ensg = var["source_feature_id"].tolist() if "source_feature_id" in var.columns else None
    aligned = align_axis_to_canonical_targets(var["gene_symbol"].astype(str).tolist(), targets, ensg)
    if not aligned["present"].all():
        missing = aligned.loc[~aligned["present"], "canonical_symbol"].head(10).tolist()
        raise ValueError(f"Missing {int((~aligned['present']).sum())} target genes in {source_label}: {missing}")
    col_indices = aligned["source_index"].to_numpy(dtype=np.int64)
    return col_indices, pd.Index(aligned["canonical_symbol"].astype(str))


def _clean_bag_obs(
    obs: pd.DataFrame,
    rows: np.ndarray,
    *,
    require_resolved_subtype: bool = True,
    source_label: str = "",
) -> pd.DataFrame:
    base_cols = [
        "donor_uid",
        "lesion_uid",
        "sample_uid",
        "tissue",
        "lineage_m2",
        "subtype_m2",
        "n_counts_visible",
        "n_genes_visible",
    ]
    extra_cols = ["myeloid_fine"] if "myeloid_fine" in obs.columns else []
    frame = obs.iloc[rows][base_cols + extra_cols].copy()
    frame.index = obs.index[rows].astype(str)

    for col in ("donor_uid", "lesion_uid", "sample_uid", "tissue", "lineage_m2", "subtype_m2", *extra_cols):
        frame[col] = frame[col].astype(str)

    if "myeloid_fine" in frame.columns:
        mf = frame["myeloid_fine"]
        fine = np.where(
            frame["lineage_m2"].eq("Mono/Mac/DC") & mf.ne("") & (~mf.isin(["nan", "None"])),
            mf,
            frame["subtype_m2"],
        )
    else:
        fine = frame["subtype_m2"].to_numpy(dtype=str)
    frame["fine_subtype"] = pd.Series(fine, index=frame.index, dtype=str)
    frame["fine_cell_type"] = frame["lineage_m2"] + "::" + frame["fine_subtype"]

    if require_resolved_subtype:
        bad_sub = set(frame["subtype_m2"]) & INVALID_SUBTYPE_LABELS
        bad_fine = set(frame["fine_subtype"]) & INVALID_SUBTYPE_LABELS
        if bad_sub or bad_fine:
            raise ValueError(
                f"Unresolved or invalid subtype labels in {source_label}: subtype_m2={sorted(bad_sub)}, fine_subtype={sorted(bad_fine)}"
            )
    return frame


def _validate_tumor_h5ad(obj, path: Path) -> tuple[pd.DataFrame, np.ndarray]:
    obs = obj.obs
    tissues = obs["tissue"].astype(str).unique().tolist()
    if len(tissues) != 1 or tissues[0] not in TUMOR_TISSUES:
        raise ValueError(f"Expected a single PT/LM tumor library at {path}, found {tissues}")
    if not obs.index.is_unique:
        raise ValueError(f"Non-unique cell identifiers in {path}")
    retained = obs["qc_retained"].to_numpy(dtype=bool)
    return obs, retained


def sample_bags(
    path: Path | str,
    *,
    min_cells: int = 50,
    max_cells: int | None = None,
    repeat_index: int = 0,
    seed: int = SEED,
    use_canonical_genes: bool = True,
    gene_universe: str | Sequence[str] | None = None,
    require_resolved_subtype: bool = True,
) -> dict[str, Bag]:
    """Load per-lineage QC-retained tumor bags (`PT`/`LM` only).

    Defaults to full bags (`max_cells=None`, `min_cells=50`). Pass `max_cells=256`
    (and optional `repeat_index`) for capped-bag sensitivity or legacy reproduction.
    """
    path = Path(path)
    obj = ad.read_h5ad(path, backed="r")
    try:
        obs, retained = _validate_tumor_h5ad(obj, path)
        col_indices, genes = _resolve_gene_slice(
            obj.var,
            use_canonical_genes=use_canonical_genes,
            gene_universe=gene_universe,
            source_label=str(path),
        )
        lin_arr = obs["lineage_m2"].astype(str).to_numpy()
        groups: dict[str, Bag] = {}
        for lineage in LINEAGES:
            cand = np.flatnonzero(retained & (lin_arr == lineage))
            population = int(len(cand))
            if population < min_cells:
                continue
            rows = select_bag_rows(
                obs.index,
                cand,
                f"{path.stem}|{lineage}",
                min_cells=min_cells,
                max_cells=max_cells,
                repeat_index=repeat_index,
                seed=seed,
            )
            frame = _clean_bag_obs(
                obs,
                rows,
                require_resolved_subtype=require_resolved_subtype,
                source_label=f"{path.name}:{lineage}",
            )
            mat = obj.X[rows, :].tocsr()
            if col_indices is not None:
                mat = mat[:, col_indices].tocsr()
            mat.eliminate_zeros()
            groups[lineage] = Bag(frame, mat, genes, population)
        return groups
    finally:
        obj.file.close()


def sample_joint_bag(
    path: Path | str,
    *,
    min_cells: int = 50,
    max_cells: int | None = None,
    repeat_index: int = 0,
    seed: int = SEED,
    lineages: Sequence[str] = LINEAGES,
    use_canonical_genes: bool = True,
    gene_universe: str | Sequence[str] | None = None,
    require_resolved_subtype: bool = True,
) -> Bag | None:
    """Load a whole-specimen multi-lineage tumor bag (`All3Lineages`, `PT`/`LM` only).

    Retains all `qc_retained` cells across `lineages` (`T/NK`, `Mono/Mac/DC`, `B/Plasma`)
    when total count `>= min_cells`, populating `obs['fine_cell_type']` for hierarchical
    multi-cell-type MIL (`HA`/`CTA`) and `MixMIL`.
    """
    path = Path(path)
    obj = ad.read_h5ad(path, backed="r")
    try:
        obs, retained = _validate_tumor_h5ad(obj, path)
        col_indices, genes = _resolve_gene_slice(
            obj.var,
            use_canonical_genes=use_canonical_genes,
            gene_universe=gene_universe,
            source_label=str(path),
        )
        lin_set = set(str(x) for x in lineages)
        lin_arr = obs["lineage_m2"].astype(str).to_numpy()
        cand = np.flatnonzero(retained & np.isin(lin_arr, list(lin_set)))
        population = int(len(cand))
        if population < min_cells:
            return None
        rows = select_bag_rows(
            obs.index,
            cand,
            f"{path.stem}|{JOINT_LINEAGE}",
            min_cells=min_cells,
            max_cells=max_cells,
            repeat_index=repeat_index,
            seed=seed,
        )
        frame = _clean_bag_obs(
            obs,
            rows,
            require_resolved_subtype=require_resolved_subtype,
            source_label=f"{path.name}:{JOINT_LINEAGE}",
        )
        mat = obj.X[rows, :].tocsr()
        if col_indices is not None:
            mat = mat[:, col_indices].tocsr()
        mat.eliminate_zeros()
        return Bag(frame, mat, genes, population)
    finally:
        obj.file.close()


def specimen_immune_composition(path: Path | str) -> pd.DataFrame:
    """Compute exact full-specimen QC-retained immune subtype proportions without subsampling noise.

    Enforces `PT`/`LM` tumor check so `NC`/`NL` normal tissues never enter tumor composition tables.
    """
    path = Path(path)
    obj = ad.read_h5ad(path, backed="r")
    try:
        obs, retained = _validate_tumor_h5ad(obj, path)
        ret_rows = np.flatnonzero(retained)
        if len(ret_rows) == 0:
            raise ValueError(f"Zero QC-retained cells in {path}")

        sub = obs.iloc[ret_rows].copy()
        sample_uid = str(sub["sample_uid"].astype(str).iloc[0])
        donor_uid = str(sub["donor_uid"].astype(str).iloc[0])
        lesion_uid = str(sub["lesion_uid"].astype(str).iloc[0])
        tissue = str(sub["tissue"].astype(str).iloc[0])

        lin_s = sub["lineage_m2"].astype(str)
        sub_s = sub["subtype_m2"].astype(str)
        if "myeloid_fine" in sub.columns:
            mf_s = sub["myeloid_fine"].astype(str)
            fine_s = pd.Series(
                np.where(
                    lin_s.eq("Mono/Mac/DC") & mf_s.ne("") & (~mf_s.isin(["nan", "None"])),
                    mf_s,
                    sub_s,
                ),
                index=sub.index,
                dtype=str,
            )
        else:
            fine_s = sub_s

        n_total_immune = int(lin_s.isin(RESOLVED_IMMUNE_LINEAGES).sum())
        n_core_three = int(lin_s.isin(LINEAGES).sum())

        core_mask = lin_s.isin(LINEAGES)
        core_df = pd.DataFrame(
            {
                "lineage_m2": lin_s[core_mask].to_numpy(dtype=str),
                "subtype_m2": sub_s[core_mask].to_numpy(dtype=str),
                "fine_subtype": fine_s[core_mask].to_numpy(dtype=str),
            }
        )
        lin_totals = core_df.groupby("lineage_m2").size().to_dict()
        grouped = (
            core_df.groupby(["lineage_m2", "fine_subtype"], as_index=False)
            .agg(subtype_m2=("subtype_m2", "first"), n_cells=("subtype_m2", "size"))
            .sort_values(["lineage_m2", "fine_subtype"])
            .reset_index(drop=True)
        )
        grouped.insert(0, "tissue", tissue)
        grouped.insert(0, "lesion_uid", lesion_uid)
        grouped.insert(0, "donor_uid", donor_uid)
        grouped.insert(0, "sample_uid", sample_uid)
        grouped["n_lineage_cells"] = grouped["lineage_m2"].map(lin_totals).fillna(0).astype(int)
        grouped["n_core_three_immune"] = n_core_three
        grouped["n_total_resolved_immune"] = n_total_immune
        grouped["prop_within_lineage"] = np.where(
            grouped["n_lineage_cells"] > 0,
            grouped["n_cells"] / grouped["n_lineage_cells"],
            0.0,
        )
        grouped["prop_core_three_immune"] = (
            grouped["n_cells"] / n_core_three if n_core_three > 0 else 0.0
        )
        grouped["prop_total_resolved_immune"] = (
            grouped["n_cells"] / n_total_immune if n_total_immune > 0 else 0.0
        )
        return grouped
    finally:
        obj.file.close()
