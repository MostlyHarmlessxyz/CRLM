"""Single-cell gene program scoring: Mann-Whitney UCell, detection-matched debiasing, and legacy rank scores.

Provides:
1. ``ucell_scores``: Exact full-axis Mann-Whitney UCell (Andreatta & Carmona 2021)
   handling both positive-expression descending midranks and measured-zero full-axis
   midranks ``(K_c + 1 + G) / 2`` in pure sparse arithmetic without densifying ``matrix``.
2. ``compute_reference_gene_stats`` & ``build_detection_matched_controls``:
   Deterministic ``n_bins=20`` expression/detection stratification on development
   reference data (GSE245552) with symmetric intra-bin caliper matching to lock
   ``n_ctrl=10`` non-member control genes per target gene for frozen cross-cohort
   reuse (never adaptively refit on external cohorts).
3. ``debiased_ucell_scores``: Depth-debiased single-cell UCell ``U'(P) - U'(P_ctrl)``,
   eliminating ``K_c < max_rank`` rank-budget inflation and discrete UMI tie cliffs.
4. ``log1p_cp10k_scores``: Scale-checked ``log1p(CP10k)`` mean program scores (raw or
   matched-control-subtracted) on a supplied common nuclear gene axis.
5. ``legacy_positive_rank_scores`` & backward-compatible ``rank_scores(..., mode=...)``:
   Bit-exact preservation of the V1 ``positive_midrank_1500`` formula alongside
   ``mode="ucell_debiased"`` and ``mode="ucell_raw"``.
"""

from __future__ import annotations

import hashlib
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

SUPPORTED_MODES = (
    "legacy_positive",
    "ucell_raw",
    "ucell_debiased",
    "log1p_cp10k_raw",
    "log1p_cp10k_debiased",
)


def _program_seed(seed: int, *tokens: object) -> int:
    payload = "|".join([str(int(seed)), *(str(t) for t in tokens)]).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _validate_axis_and_programs(
    n_cols: int,
    symbols: Sequence[str],
    programs: Mapping[str, Sequence[str]],
    *,
    max_rank: int | None = None,
    require_size_le_max_rank: bool = False,
) -> tuple[pd.Index, list[str], sparse.csr_matrix, np.ndarray]:
    sym_index = pd.Index([str(s) for s in symbols])
    if len(sym_index) == 0 or not sym_index.is_unique or len(sym_index) != n_cols:
        raise ValueError("Require a nonempty unique gene axis matching matrix columns")
    if max_rank is not None and (int(max_rank) != max_rank or max_rank < 1):
        raise ValueError("Require a positive integer max_rank")
    names = list(programs)
    genes_in_sets: list[int] = []
    set_ids: list[int] = []
    weights: list[float] = []
    sizes = np.empty(len(names), dtype=np.float64)
    for number, name in enumerate(names):
        members = [str(g) for g in programs[name]]
        indices = sym_index.get_indexer(members)
        if not members or len(set(members)) != len(members) or (indices < 0).any():
            raise ValueError(f"{name}: supply unique, nonempty measured program members")
        if require_size_le_max_rank and max_rank is not None and len(members) > max_rank:
            raise ValueError(
                f"{name}: signature length ({len(members)}) exceeds max_rank ({max_rank}) for Mann-Whitney UCell"
            )
        sizes[number] = float(len(members))
        genes_in_sets.extend(indices.tolist())
        set_ids.extend([number] * len(indices))
        weights.extend([1.0 / len(indices)] * len(indices))
    design = sparse.csr_matrix(
        (weights, (genes_in_sets, set_ids)),
        shape=(len(sym_index), len(names)),
        dtype=np.float64,
    )
    return sym_index, names, design, sizes


def legacy_positive_rank_scores(
    matrix,
    symbols: Sequence[str],
    programs: Mapping[str, Sequence[str]],
    *,
    max_rank: int = 1500,
    chunk_size: int = 512,
) -> pd.DataFrame:
    """Bit-exact V1 legacy positive-expression descending midrank score (`positive_midrank_1500`).

    Positive values receive `max(0, 1 - (descending_midrank - 1) / max_rank)` and
    measured zeros receive zero weight. Preserved strictly for Track 1 (`techfix`)
    backward reproduction; use `ucell_scores` or `debiased_ucell_scores` for V1.1.
    """
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    _, names, design, _ = _validate_axis_and_programs(
        matrix.shape[1], symbols, programs, max_rank=max_rank, require_size_le_max_rank=False
    )
    if not names:
        return pd.DataFrame(index=np.arange(matrix.shape[0]), columns=names, dtype=np.float64)
    relevant = np.diff(design.indptr) > 0
    scores = np.zeros((matrix.shape[0], len(names)), dtype=np.float64)
    for start in range(0, matrix.shape[0], chunk_size):
        chunk = sparse.csr_matrix(matrix[start : start + chunk_size], copy=True)
        chunk.sum_duplicates()
        chunk.eliminate_zeros()
        if not np.isfinite(chunk.data).all() or (chunk.data < 0).any():
            raise ValueError("Rank inputs must have finite nonnegative expression")
        ranks = np.zeros(chunk.nnz, dtype=np.float64)
        for a, b in zip(chunk.indptr[:-1], chunk.indptr[1:]):
            values = chunk.data[a:b]
            keep = relevant[chunk.indices[a:b]]
            if not keep.any():
                continue
            ordered = np.sort(values)
            target = values[keep]
            left = np.searchsorted(ordered, target, side="left")
            right = np.searchsorted(ordered, target, side="right")
            descending_midrank = len(ordered) + 0.5 - 0.5 * (left + right)
            ranks[a:b][keep] = np.maximum(0.0, 1.0 - (descending_midrank - 1.0) / max_rank)
        chunk.data = ranks
        chunk.eliminate_zeros()
        scores[start : start + chunk.shape[0]] = (chunk @ design).toarray()
    return pd.DataFrame(scores, columns=names)


def ucell_scores(
    matrix,
    symbols: Sequence[str],
    programs: Mapping[str, Sequence[str]],
    *,
    max_rank: int = 1500,
    chunk_size: int = 512,
) -> pd.DataFrame:
    """Exact full-axis Mann-Whitney UCell score (Andreatta & Carmona 2021).

    Across all `G = len(symbols)` measured genes in cell `c` with `K_c` positive
    genes:
    - Positive genes receive descending midranks `r_{c,g} in [1, K_c]`.
    - Measured zero genes receive their exact tied full-axis midrank
      `r_zero(c) = (K_c + 1 + G) / 2` (exact even when `G < max_rank`).
    - All ranks `r > max_rank` are capped to `r_tilde = max_rank + 1`.
    - Returns `U'_c(S) = max(0, 1 - (sum_{g in S} r_tilde_{c,g} - m*(m+1)/2) / (m * max_rank))`
      with unit upper bound `1.0` independent of signature size `m`.
    """
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    sym_index, names, design, sizes = _validate_axis_and_programs(
        matrix.shape[1], symbols, programs, max_rank=max_rank, require_size_le_max_rank=True
    )
    n_rows = matrix.shape[0]
    if not names:
        return pd.DataFrame(index=np.arange(n_rows), columns=names, dtype=np.float64)
    g_total = float(len(sym_index))
    cap_val = float(max_rank) + 1.0
    self_midrank = 0.5 * (sizes + 1.0)
    relevant = np.diff(design.indptr) > 0
    scores = np.zeros((n_rows, len(names)), dtype=np.float64)

    for start in range(0, n_rows, chunk_size):
        chunk = sparse.csr_matrix(matrix[start : start + chunk_size], copy=True, dtype=np.float64)
        chunk.sum_duplicates()
        chunk.eliminate_zeros()
        if not np.isfinite(chunk.data).all() or (chunk.data < 0).any():
            raise ValueError("UCell inputs must have finite nonnegative expression")
        n_sub = chunk.shape[0]
        deltas = np.zeros(chunk.nnz, dtype=np.float64)
        r0_vec = np.empty(n_sub, dtype=np.float64)
        for row_idx, (a, b) in enumerate(zip(chunk.indptr[:-1], chunk.indptr[1:])):
            k_pos = b - a
            r_zero = 0.5 * (k_pos + 1.0 + g_total)
            r0 = r_zero if r_zero <= max_rank else cap_val
            r0_vec[row_idx] = r0
            if k_pos == 0:
                continue
            keep = relevant[chunk.indices[a:b]]
            if not keep.any():
                continue
            values = chunk.data[a:b]
            ordered = np.sort(values)
            target = values[keep]
            left = np.searchsorted(ordered, target, side="left")
            right = np.searchsorted(ordered, target, side="right")
            descending_midrank = k_pos + 0.5 - 0.5 * (left + right)
            capped = np.where(descending_midrank <= max_rank, descending_midrank, cap_val)
            deltas[a:b][keep] = r0 - capped
        chunk.data = deltas
        chunk.eliminate_zeros()
        mean_capped_rank = r0_vec[:, None] - (chunk @ design).toarray()
        raw_ucell = 1.0 - (mean_capped_rank - self_midrank[None, :]) / float(max_rank)
        scores[start : start + n_sub] = np.clip(raw_ucell, 0.0, 1.0)

    return pd.DataFrame(scores, columns=names)


def _iter_matrices(matrices) -> Iterable:
    if sparse.issparse(matrices) or isinstance(matrices, (np.ndarray, pd.DataFrame)):
        yield matrices
    else:
        yield from matrices


def compute_reference_gene_stats(
    matrices,
    symbols: Sequence[str],
    *,
    input_scale: str = "counts",
    scale_factor: float = 10000.0,
    n_bins: int = 20,
    chunk_size: int = 1024,
) -> pd.DataFrame:
    """Compute per-gene detection rate, mean CP10k, mean log1p(CP10k), and expression bins.

    Intended to be called ONLY on development cohort (`GSE245552`) cells so that
    bin assignments and matched control sets can be frozen prior to external validation.
    Accepts either a single `(cells x genes)` matrix or an iterable of matrices
    aligned to `symbols`.
    """
    if input_scale not in ("counts", "log1p"):
        raise ValueError("input_scale must be 'counts' or 'log1p'")
    if scale_factor <= 0 or not np.isfinite(scale_factor):
        raise ValueError("scale_factor must be finite and positive")
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    sym_index = pd.Index([str(s) for s in symbols])
    if len(sym_index) == 0 or not sym_index.is_unique:
        raise ValueError("Require a nonempty unique gene axis")

    n_genes = len(sym_index)
    detected_counts = np.zeros(n_genes, dtype=np.int64)
    sum_cp10k = np.zeros(n_genes, dtype=np.float64)
    sum_log1p = np.zeros(n_genes, dtype=np.float64)
    total_cells = 0

    for mat in _iter_matrices(matrices):
        if mat.shape[1] != n_genes:
            raise ValueError(f"Matrix columns ({mat.shape[1]}) do not match symbols ({n_genes})")
        for start in range(0, mat.shape[0], chunk_size):
            chunk = sparse.csr_matrix(mat[start : start + chunk_size], copy=True, dtype=np.float64)
            chunk.sum_duplicates()
            chunk.eliminate_zeros()
            if not np.isfinite(chunk.data).all() or (chunk.data < 0).any():
                raise ValueError("Reference expression must be finite and nonnegative")
            if input_scale == "log1p":
                chunk.data = np.expm1(chunk.data)
            detected_counts += np.bincount(chunk.indices, minlength=n_genes)
            row_sums = np.asarray(chunk.sum(axis=1)).ravel()
            inv = np.where(row_sums > 0.0, scale_factor / row_sums, 0.0)
            scaled = chunk.multiply(inv[:, None]).tocsr()
            sum_cp10k += np.asarray(scaled.sum(axis=0)).ravel()
            scaled.data = np.log1p(scaled.data)
            sum_log1p += np.asarray(scaled.sum(axis=0)).ravel()
            total_cells += chunk.shape[0]

    if total_cells == 0:
        raise ValueError("Reference matrices contain zero cells")

    stats = pd.DataFrame(
        {
            "gene_symbol": sym_index.to_numpy(dtype=str),
            "n_cells_reference": int(total_cells),
            "n_detected": detected_counts,
            "detection_rate": detected_counts.astype(np.float64) / float(total_cells),
            "mean_cp10k": sum_cp10k / float(total_cells),
            "mean_log1p_cp10k": sum_log1p / float(total_cells),
        }
    )
    stats["expression_bin"] = assign_expression_bins(stats, n_bins=n_bins).to_numpy()
    return stats


def _prepare_stats_frame(reference_gene_stats: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(reference_gene_stats, pd.DataFrame):
        raise TypeError("reference_gene_stats must be a pandas DataFrame")
    frame = reference_gene_stats.copy()
    if "gene_symbol" not in frame.columns:
        if frame.index.name == "gene_symbol" or frame.index.dtype == object:
            frame = frame.reset_index().rename(columns={frame.index.name or "index": "gene_symbol"})
        else:
            raise ValueError("reference_gene_stats must include a 'gene_symbol' column or string index")
    frame["gene_symbol"] = frame["gene_symbol"].astype(str)
    if frame["gene_symbol"].duplicated().any():
        raise ValueError("reference_gene_stats has duplicate gene_symbol entries")
    if "detection_rate" not in frame.columns:
        raise ValueError("reference_gene_stats requires a 'detection_rate' column")
    if "mean_log1p_cp10k" not in frame.columns:
        frame["mean_log1p_cp10k"] = frame.get("mean_cp10k", frame["detection_rate"]).astype(np.float64)
    return frame.set_index("gene_symbol", drop=False)


def assign_expression_bins(
    reference_gene_stats: pd.DataFrame,
    *,
    candidate_universe: Sequence[str] | None = None,
    n_bins: int = 20,
) -> pd.Series:
    """Assign deterministic expression/detection bins (`0 .. n_bins - 1`).

    Unexpressed genes (`detection_rate == 0`) occupy bin 0 when present; expressed
    genes are stratified by `(detection_rate, mean_log1p_cp10k, gene_symbol)` using
    tail-refined quantiles `1 - sqrt(1 - u)` so both sparse and highly expressed
    transcriptome regimes have tight intra-bin detection widths.
    """
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    frame = _prepare_stats_frame(reference_gene_stats)
    if candidate_universe is not None:
        universe = pd.Index([str(s) for s in candidate_universe])
        if not universe.is_unique or not universe.isin(frame.index).all():
            missing = sorted(set(universe) - set(frame.index))[:10]
            raise ValueError(f"candidate_universe must be unique and covered by reference_gene_stats; missing={missing}")
        frame = frame.loc[universe]

    syms = frame["gene_symbol"].to_numpy(dtype=str)
    det = frame["detection_rate"].to_numpy(dtype=np.float64)
    mlp = frame["mean_log1p_cp10k"].to_numpy(dtype=np.float64)
    if not np.isfinite(det).all() or not np.isfinite(mlp).all() or (det < 0).any():
        raise ValueError("reference_gene_stats contains nonfinite or negative values")

    bins = np.zeros(len(frame), dtype=np.int32)
    if n_bins == 1 or len(frame) == 0:
        return pd.Series(bins, index=frame.index, name="expression_bin")

    zero_mask = (det <= 0.0) & (mlp <= 0.0)
    if zero_mask.any() and (~zero_mask).any():
        pos_idx = np.flatnonzero(~zero_mask)
        order = pos_idx[np.lexsort((syms[pos_idx], mlp[pos_idx], det[pos_idx]))]
        u = np.arange(len(order), dtype=np.float64) / float(len(order))
        n_pos_bins = n_bins - 1
        bins[order] = 1 + np.minimum(
            n_pos_bins - 1, np.floor(n_pos_bins * (1.0 - np.sqrt(1.0 - u))).astype(np.int32)
        )
    else:
        order = np.lexsort((syms, mlp, det))
        u = np.arange(len(order), dtype=np.float64) / float(len(order))
        bins[order] = np.minimum(
            n_bins - 1, np.floor(n_bins * (1.0 - np.sqrt(1.0 - u))).astype(np.int32)
        )
    return pd.Series(bins, index=frame.index, name="expression_bin")


def build_detection_matched_controls(
    programs: Mapping[str, Sequence[str]],
    reference_gene_stats: pd.DataFrame,
    *,
    candidate_universe: Sequence[str] | None = None,
    exclude_all_program_genes: bool = False,
    n_bins: int = 20,
    n_ctrl: int = 10,
    seed: int = 20261121,
) -> dict[str, dict[str, list[str]]]:
    """Deterministically select `n_ctrl` non-member same-bin background genes per target gene.

    Returns `{program_name: {target_gene: [ctrl_gene_0, ..., ctrl_gene_{n_ctrl-1}]}}`.
    Selection is seeded independently per `(seed, program_name)` so adding/removing
    other programs never perturbs a program's frozen control set. Within each bin,
    non-member genes are ordered by `(detection_rate, mean_log1p_cp10k, gene_symbol)`
    and sampled without replacement from a centered symmetric caliper window while
    preventing duplicate control genes within each of the `n_ctrl` control replicates.
    """
    if n_bins < 1 or n_ctrl < 1:
        raise ValueError("n_bins and n_ctrl must be >= 1")
    frame = _prepare_stats_frame(reference_gene_stats)
    universe = (
        pd.Index([str(s) for s in candidate_universe])
        if candidate_universe is not None
        else pd.Index(frame.index)
    )
    if len(universe) == 0 or not universe.is_unique or not universe.isin(frame.index).all():
        raise ValueError("candidate_universe must be nonempty, unique, and present in reference_gene_stats")
    sub = frame.loc[universe].copy()
    sub["bin"] = assign_expression_bins(sub, n_bins=n_bins).to_numpy()

    syms = sub["gene_symbol"].to_numpy(dtype=str)
    det = sub["detection_rate"].to_numpy(dtype=np.float64)
    mlp = sub["mean_log1p_cp10k"].to_numpy(dtype=np.float64)
    bins = sub["bin"].to_numpy(dtype=np.int32)
    sym_to_pos = {s: i for i, s in enumerate(syms)}

    global_excluded = (
        {str(g) for members in programs.values() for g in members}
        if exclude_all_program_genes
        else set()
    )

    controls: dict[str, dict[str, list[str]]] = {}
    for pname, raw_members in programs.items():
        members = [str(g) for g in raw_members]
        if not members or len(set(members)) != len(members):
            raise ValueError(f"{pname}: program members must be nonempty and unique")
        missing = [g for g in members if g not in sym_to_pos]
        if missing:
            raise ValueError(f"{pname}: program members absent from reference universe: {missing[:10]}")

        forbidden = set(members) | global_excluded
        non_member_mask = np.array([s not in forbidden for s in syms], dtype=bool)
        if not non_member_mask.any():
            raise ValueError(f"{pname}: zero non-member background genes available in universe")

        rng = np.random.default_rng(_program_seed(seed, pname))
        used_in_rep: list[set[str]] = [set() for _ in range(n_ctrl)]
        prog_map: dict[str, list[str]] = {}

        # Sort target genes deterministically so replicate uniqueness is order-invariant
        ordered_members = sorted(
            members,
            key=lambda g: (bins[sym_to_pos[g]], det[sym_to_pos[g]], mlp[sym_to_pos[g]], g),
        )
        temp_map: dict[str, list[str]] = {}
        half_win = max(n_ctrl + 5, 15)

        for g in ordered_members:
            gi = sym_to_pos[g]
            b_id = bins[gi]
            pool_idx = np.flatnonzero((bins == b_id) & non_member_mask)
            if len(pool_idx) == 0:
                pool_idx = np.flatnonzero(non_member_mask)
            pool_ord = pool_idx[
                np.lexsort((syms[pool_idx], mlp[pool_idx], det[pool_idx]))
            ]
            pos = int(np.searchsorted(det[pool_ord], det[gi], side="left"))
            # Centered symmetric caliper around `pos` inside the same bin
            sym_radius = min(half_win, pos, len(pool_ord) - pos)
            if 2 * sym_radius >= max(2 * n_ctrl, 12):
                win_indices = pool_ord[pos - sym_radius : pos + sym_radius]
            elif len(pool_ord) > 2 * half_win:
                start = min(max(0, pos - half_win), len(pool_ord) - 2 * half_win)
                win_indices = pool_ord[start : start + 2 * half_win]
            else:
                win_indices = pool_ord

            win_syms = syms[win_indices]
            all_bin_syms = syms[pool_ord]
            chosen_for_g: list[str] = []
            used_for_this_gene: set[str] = set()

            # Permute candidate lists once per target gene using deterministic RNG
            perm_win = rng.permutation(win_syms).tolist()
            perm_bin = rng.permutation(all_bin_syms).tolist()

            for rep_j in range(n_ctrl):
                pick = next(
                    (c for c in perm_win if c not in used_for_this_gene and c not in used_in_rep[rep_j]),
                    None,
                )
                if pick is None:
                    pick = next(
                        (c for c in perm_bin if c not in used_for_this_gene and c not in used_in_rep[rep_j]),
                        None,
                    )
                if pick is None:
                    pick = next((c for c in perm_bin if c not in used_for_this_gene), None)
                if pick is None:
                    pick = perm_bin[rep_j % len(perm_bin)]
                chosen_for_g.append(str(pick))
                used_for_this_gene.add(str(pick))
                used_in_rep[rep_j].add(str(pick))

            temp_map[g] = chosen_for_g

        # Preserve caller's member order in returned dictionary
        for g in members:
            prog_map[g] = temp_map[g]
        controls[pname] = prog_map

    return controls


def controls_to_dataframe(
    matched_controls: Mapping[str, Mapping[str, Sequence[str]]],
    reference_gene_stats: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Serialize `build_detection_matched_controls` output to a flat tidy DataFrame."""
    stats = _prepare_stats_frame(reference_gene_stats) if reference_gene_stats is not None else None
    rows = []
    for pname, gmap in matched_controls.items():
        for target_gene, ctrls in gmap.items():
            for rep_idx, ctrl_gene in enumerate(ctrls):
                entry = {
                    "program": str(pname),
                    "target_gene": str(target_gene),
                    "replicate": int(rep_idx),
                    "control_gene": str(ctrl_gene),
                }
                if stats is not None and target_gene in stats.index and ctrl_gene in stats.index:
                    entry["target_detection_rate"] = float(stats.at[target_gene, "detection_rate"])
                    entry["control_detection_rate"] = float(stats.at[ctrl_gene, "detection_rate"])
                    if "expression_bin" in stats.columns:
                        entry["expression_bin"] = int(stats.at[target_gene, "expression_bin"])
                rows.append(entry)
    return pd.DataFrame(rows)


def controls_from_dataframe(frame: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    """Reconstruct frozen `matched_controls` dictionary from a saved DataFrame."""
    required = {"program", "target_gene", "replicate", "control_gene"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Control table missing columns: {sorted(missing)}")
    ordered = frame.sort_values(["program", "replicate"], kind="mergesort")
    out: dict[str, dict[str, list[str]]] = {}
    for pname, pgrp in ordered.groupby("program", sort=False):
        gmap: dict[str, list[str]] = {}
        for tgene, tgrp in pgrp.sort_values("replicate", kind="mergesort").groupby("target_gene", sort=False):
            gmap[str(tgene)] = [str(c) for c in tgrp["control_gene"].tolist()]
        out[str(pname)] = gmap
    return out


def _resolve_control_replicates(
    programs: Mapping[str, Sequence[str]],
    symbols_set: set[str],
    matched_controls: Mapping,
    *,
    min_ctrl_fraction: float = 0.8,
) -> dict[str, list[list[str]]]:
    """Normalize supported `matched_controls` structures into measured replicate lists."""
    replicates: dict[str, list[list[str]]] = {}
    for pname, raw_members in programs.items():
        if pname not in matched_controls:
            raise ValueError(f"{pname}: missing from supplied matched_controls")
        members = [str(g) for g in raw_members]
        spec = matched_controls[pname]

        if isinstance(spec, Mapping):
            missing_targets = [g for g in members if g not in spec]
            if missing_targets:
                raise ValueError(f"{pname}: matched_controls missing target genes {missing_targets[:10]}")
            n_rep = max(len(spec[g]) for g in members)
            if n_rep < 1:
                raise ValueError(f"{pname}: empty control list in matched_controls")
            reps: list[list[str]] = []
            total_expected = len(members) * n_rep
            total_measured = 0
            for j in range(n_rep):
                seen: set[str] = set()
                rep_genes: list[str] = []
                for g in members:
                    g_ctrls = spec[g]
                    if j >= len(g_ctrls):
                        continue
                    cg = str(g_ctrls[j])
                    if cg in symbols_set:
                        total_measured += 1
                        if cg not in seen:
                            seen.add(cg)
                            rep_genes.append(cg)
                if not rep_genes:
                    raise ValueError(f"{pname}: control replicate {j} has zero measured genes on axis")
                reps.append(rep_genes)
            if total_measured < total_expected * min_ctrl_fraction:
                raise ValueError(
                    f"{pname}: insufficient measured frozen control genes ({total_measured}/{total_expected})"
                )
            replicates[pname] = reps
        elif isinstance(spec, Sequence) and not isinstance(spec, (str, bytes)):
            items = list(spec)
            if not items:
                raise ValueError(f"{pname}: empty matched_controls entry")
            if isinstance(items[0], Sequence) and not isinstance(items[0], (str, bytes)):
                reps = []
                total_expected = sum(len(r) for r in items)
                total_measured = 0
                for j, raw_rep in enumerate(items):
                    seen = set()
                    rep_genes = []
                    for cg in raw_rep:
                        s_cg = str(cg)
                        if s_cg in symbols_set:
                            total_measured += 1
                            if s_cg not in seen:
                                seen.add(s_cg)
                                rep_genes.append(s_cg)
                    if not rep_genes:
                        raise ValueError(f"{pname}: control replicate {j} has zero measured genes on axis")
                    reps.append(rep_genes)
                if total_expected > 0 and total_measured < total_expected * min_ctrl_fraction:
                    raise ValueError(
                        f"{pname}: insufficient measured frozen control genes ({total_measured}/{total_expected})"
                    )
                replicates[pname] = reps
            else:
                flat = [str(cg) for cg in items]
                measured = list(dict.fromkeys(cg for cg in flat if cg in symbols_set))
                if not measured or len(measured) < len(set(flat)) * min_ctrl_fraction:
                    raise ValueError(f"{pname}: insufficient measured frozen control genes")
                m_target = len(members)
                n_chunks = max(1, int(np.ceil(len(measured) / m_target)))
                replicates[pname] = [list(arr) for arr in np.array_split(measured, n_chunks) if len(arr)]
        else:
            raise TypeError(f"{pname}: unsupported matched_controls value type {type(spec)!r}")
    return replicates


def debiased_ucell_scores(
    matrix,
    symbols: Sequence[str],
    programs: Mapping[str, Sequence[str]],
    *,
    matched_controls=None,
    reference_gene_stats: pd.DataFrame | None = None,
    max_rank: int = 1500,
    n_bins: int = 20,
    n_ctrl: int = 10,
    seed: int = 20261121,
    min_ctrl_fraction: float = 0.8,
    chunk_size: int = 512,
    return_components: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute depth-debiased UCell scores `U'_c(P) - U'_c(P_ctrl)` in a single sparse pass.

    Requires either pre-frozen `matched_controls` (mandatory for external cohorts)
    or development `reference_gene_stats` (for discovery cohort construction).
    Adaptive self-fitting on unlabeled/external `matrix` is strictly prohibited.
    """
    if matched_controls is None and reference_gene_stats is None:
        raise ValueError(
            "debiased_ucell_scores requires pre-frozen matched_controls or development "
            "reference_gene_stats; adaptive self-fitting on the target matrix is prohibited"
        )
    sym_list = [str(s) for s in symbols]
    sym_set = set(sym_list)
    if len(sym_set) != len(sym_list) or len(sym_list) != matrix.shape[1]:
        raise ValueError("Require a nonempty unique gene axis matching matrix columns")

    if matched_controls is None:
        matched_controls = build_detection_matched_controls(
            programs,
            reference_gene_stats,
            candidate_universe=sym_list,
            n_bins=n_bins,
            n_ctrl=n_ctrl,
            seed=seed,
        )

    names = list(programs)
    if not names:
        empty = pd.DataFrame(index=np.arange(matrix.shape[0]), columns=names, dtype=np.float64)
        return (empty, empty.copy(), empty.copy()) if return_components else empty

    ctrl_replicates = _resolve_control_replicates(
        programs, sym_set, matched_controls, min_ctrl_fraction=min_ctrl_fraction
    )

    combined: dict[str, Sequence[str]] = {}
    ctrl_keys: dict[str, list[str]] = {}
    for pname in names:
        combined[pname] = programs[pname]
        keys_p: list[str] = []
        for j, rep_genes in enumerate(ctrl_replicates[pname]):
            key = f"__ctrl__{pname}__{j}"
            while key in combined or key in programs:
                key = f"_{key}"
            combined[key] = rep_genes
            keys_p.append(key)
        ctrl_keys[pname] = keys_p

    all_ucell = ucell_scores(matrix, sym_list, combined, max_rank=max_rank, chunk_size=chunk_size)
    raw_df = all_ucell[names].copy()
    ctrl_df = pd.DataFrame(
        {pname: all_ucell[ctrl_keys[pname]].mean(axis=1).to_numpy(dtype=np.float64) for pname in names},
        index=raw_df.index,
    )
    debiased_df = raw_df - ctrl_df
    if return_components:
        return debiased_df, raw_df, ctrl_df
    return debiased_df


def log1p_cp10k_scores(
    matrix,
    symbols: Sequence[str],
    programs: Mapping[str, Sequence[str]],
    *,
    matched_controls=None,
    reference_gene_stats: pd.DataFrame | None = None,
    debiased: bool | None = None,
    input_scale: str = "counts",
    input_is_log1p: bool = False,
    renormalize_subaxis: bool = True,
    scale_factor: float = 10000.0,
    n_bins: int = 20,
    n_ctrl: int = 10,
    seed: int = 20261121,
    min_ctrl_fraction: float = 0.8,
    chunk_size: int = 512,
    return_components: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score `log1p(CP10k)` program means (raw or matched-control-subtracted) on a common axis.

    - When `input_scale="counts"` (default), normalizes each cell's row sum on
      `symbols` to `scale_factor` (`10000.0`) and applies `log1p`.
    - When `input_scale="log1p"` (or `input_is_log1p=True`, e.g. GSE164522
      published `log1p(CP10k)`), recovers linear proportions via `expm1` on the
      supplied `symbols` axis and re-normalizes to `scale_factor` when
      `renormalize_subaxis=True` so missing non-nuclear genes do not shift scale.
    - When `matched_controls` or `reference_gene_stats` is supplied (or
      `debiased=True`), subtracts the matched-control `log1p(CP10k)` mean.
    """
    if input_is_log1p:
        input_scale = "log1p"
    if input_scale not in ("counts", "log1p"):
        raise ValueError("input_scale must be 'counts' or 'log1p'")
    if scale_factor <= 0 or not np.isfinite(scale_factor):
        raise ValueError("scale_factor must be finite and positive")
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")

    use_debiased = (
        (matched_controls is not None or reference_gene_stats is not None)
        if debiased is None
        else bool(debiased)
    )
    if use_debiased and matched_controls is None and reference_gene_stats is None:
        raise ValueError(
            "Debiased log1p_cp10k_scores requires frozen matched_controls or development "
            "reference_gene_stats; adaptive self-fitting on target matrices is prohibited"
        )

    sym_list = [str(s) for s in symbols]
    sym_set = set(sym_list)
    names = list(programs)
    if not names:
        empty = pd.DataFrame(index=np.arange(matrix.shape[0]), columns=names, dtype=np.float64)
        return (empty, empty.copy(), empty.copy()) if return_components else empty

    combined: dict[str, Sequence[str]] = {}
    ctrl_keys: dict[str, list[str]] = {}
    if use_debiased:
        if matched_controls is None:
            matched_controls = build_detection_matched_controls(
                programs,
                reference_gene_stats,
                candidate_universe=sym_list,
                n_bins=n_bins,
                n_ctrl=n_ctrl,
                seed=seed,
            )
        ctrl_replicates = _resolve_control_replicates(
            programs, sym_set, matched_controls, min_ctrl_fraction=min_ctrl_fraction
        )
        for pname in names:
            combined[pname] = programs[pname]
            keys_p: list[str] = []
            for j, rep_genes in enumerate(ctrl_replicates[pname]):
                key = f"__ctrl__{pname}__{j}"
                while key in combined or key in programs:
                    key = f"_{key}"
                combined[key] = rep_genes
                keys_p.append(key)
            ctrl_keys[pname] = keys_p
    else:
        for pname in names:
            combined[pname] = programs[pname]

    _, all_names, design, _ = _validate_axis_and_programs(matrix.shape[1], sym_list, combined)
    scores = np.zeros((matrix.shape[0], len(all_names)), dtype=np.float64)
    for start in range(0, matrix.shape[0], chunk_size):
        chunk = sparse.csr_matrix(matrix[start : start + chunk_size], copy=True, dtype=np.float64)
        chunk.sum_duplicates()
        chunk.eliminate_zeros()
        if not np.isfinite(chunk.data).all() or (chunk.data < 0).any():
            raise ValueError("Expression inputs must be finite and nonnegative")
        if input_scale == "counts":
            row_sums = np.asarray(chunk.sum(axis=1)).ravel()
            inv = np.where(row_sums > 0.0, scale_factor / row_sums, 0.0)
            chunk = chunk.multiply(inv[:, None]).tocsr()
            chunk.data = np.log1p(chunk.data)
        elif renormalize_subaxis:
            chunk.data = np.expm1(chunk.data)
            row_sums = np.asarray(chunk.sum(axis=1)).ravel()
            inv = np.where(row_sums > 0.0, scale_factor / row_sums, 0.0)
            chunk = chunk.multiply(inv[:, None]).tocsr()
            chunk.data = np.log1p(chunk.data)
        scores[start : start + chunk.shape[0]] = (chunk @ design).toarray()

    all_df = pd.DataFrame(scores, columns=all_names)
    raw_df = all_df[names].copy()
    if not use_debiased:
        if return_components:
            zeros = pd.DataFrame(0.0, index=raw_df.index, columns=names)
            return raw_df, raw_df.copy(), zeros
        return raw_df

    ctrl_df = pd.DataFrame(
        {pname: all_df[ctrl_keys[pname]].mean(axis=1).to_numpy(dtype=np.float64) for pname in names},
        index=raw_df.index,
    )
    debiased_df = raw_df - ctrl_df
    if return_components:
        return debiased_df, raw_df, ctrl_df
    return debiased_df


def rank_scores(
    matrix,
    symbols: Sequence[str],
    programs: Mapping[str, Sequence[str]],
    *,
    max_rank: int = 1500,
    mode: str = "legacy_positive",
    matched_controls=None,
    reference_gene_stats: pd.DataFrame | None = None,
    input_scale: str = "counts",
    n_bins: int = 20,
    n_ctrl: int = 10,
    seed: int = 20261121,
    min_ctrl_fraction: float = 0.8,
    chunk_size: int = 512,
) -> pd.DataFrame:
    """Unified program scoring entry point preserving Track 1 legacy call compatibility.

    Modes:
    - ``"legacy_positive"`` (default for legacy 3-argument callers): bit-exact V1
      `positive_midrank_1500` implementation (`legacy_positive_rank_scores`).
    - ``"ucell_raw"``: exact full-axis Mann-Whitney UCell (`ucell_scores`).
    - ``"ucell_debiased"``: detection-matched control-subtracted UCell
      (`debiased_ucell_scores`; requires `matched_controls` or `reference_gene_stats`).
    - ``"log1p_cp10k_raw"`` / ``"log1p_cp10k_debiased"``: scale-normalized
      `log1p(CP10k)` mean program score (raw or matched-control-subtracted).
    """
    if mode not in SUPPORTED_MODES:
        raise ValueError(f"Unsupported scoring mode {mode!r}; expected one of {SUPPORTED_MODES}")
    if mode == "legacy_positive":
        return legacy_positive_rank_scores(
            matrix, symbols, programs, max_rank=max_rank, chunk_size=chunk_size
        )
    if mode == "ucell_raw":
        return ucell_scores(
            matrix, symbols, programs, max_rank=max_rank, chunk_size=chunk_size
        )
    if mode == "ucell_debiased":
        return debiased_ucell_scores(
            matrix,
            symbols,
            programs,
            matched_controls=matched_controls,
            reference_gene_stats=reference_gene_stats,
            max_rank=max_rank,
            n_bins=n_bins,
            n_ctrl=n_ctrl,
            seed=seed,
            min_ctrl_fraction=min_ctrl_fraction,
            chunk_size=chunk_size,
        )
    if mode == "log1p_cp10k_raw":
        return log1p_cp10k_scores(
            matrix,
            symbols,
            programs,
            debiased=False,
            input_scale=input_scale,
            chunk_size=chunk_size,
        )
    return log1p_cp10k_scores(
        matrix,
        symbols,
        programs,
        matched_controls=matched_controls,
        reference_gene_stats=reference_gene_stats,
        debiased=True,
        input_scale=input_scale,
        n_bins=n_bins,
        n_ctrl=n_ctrl,
        seed=seed,
        min_ctrl_fraction=min_ctrl_fraction,
        chunk_size=chunk_size,
    )
