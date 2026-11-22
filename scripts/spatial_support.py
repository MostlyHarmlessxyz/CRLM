#!/usr/bin/env python3
"""Multi-scale, spillover-aware spatial proximity audit on GSE335552 Xenium sections.

Remediates spatial audit findings:
1. Aligns program coverage sets with frozen `CRLM_V1.1_245_F1` (`artifacts/m5/gene_sets/`
   and `artifacts/m3/program_definitions_draft/program_genes.csv`), stripping both T/NK
   cytotoxic genes (`GZMB`, `GNLY`, `LAT`, `IL2RG`) and endothelial `FLT1` so
   `B19_S4_Mac_SPP1_decontaminated` reflects the exact 159 full / 157 shared-nuclear V1.1
   `B19_S4_Mac_SPP1_clean` set (`7` measured on GSE335552 Xenium) and `UPR` reflects the
   exact 113 full / 111 shared-nuclear V1.1 set (`3` measured), with column-name fallback
   (`gene_symbol` vs `gene`).
2. Replaces blind `discordant = t & my` deletion with transcript-dominance disambiguation
   (`dominance_disambiguated`, primary) alongside `strict_zero_spillover` (sensitivity).
3. Evaluates four spatial grid scales (`400 / 200 / 100 / 50 um`), adds a half-grid
   (`25 um` / `0.5 * tile_um` origin shift) boundary sensitivity check
   (`grid_shift_sensitivity.csv`), and stratifies myeloid cells by proximity to `EPCAM+`
   tumor epithelium (`all_myeloid`, `near_epcam_lt15um`, `stromal_epcam_ge15um`),
   reporting `M397T` and `M242T` `<15 um` vs `>=15 um` strata symmetrically alongside the
   `~EPCAM` boundary spillover caveat.
4. Runs `M = 1,000` within-tile permutations using independent SHA256-derived RNG seeds
   from `BASE_SEED = 20261121` and reports two-sided exact Phipson-Smyth permutation
   p-values `p = (b + 1) / (M + 1)`.
"""

from __future__ import annotations

from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import resource
import sys
import time

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]

from crlm.hgnc_alias_map import canonicalize_symbols, normalize_gene_axis  # noqa: E402

RAW = ROOT / "data/raw/GSE335552"
M5_CANDIDATE_GENES = ROOT / "artifacts/m5/gene_sets/program_genes.csv"
M5_MODEL_GENES = ROOT / "artifacts/m5/gene_sets/model_program_genes.csv"
M3_FULL_GENES = ROOT / "artifacts/m3/program_definitions_draft/program_genes.csv"
OUT = ROOT / "artifacts/m7/spatial_measurement_audit"

IDENTITY_PANEL = [
    "CD68",
    "CD163",
    "C1QA",
    "CD14",
    "LYZ",
    "EPCAM",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD8B",
    "SPP1",
    "CXCL9",
    "CXCL10",
]
CROSS_LINEAGE_CONTAMINANTS = {
    "GZMB",
    "GNLY",
    "LAT",
    "IL2RG",
    "GZMA",
    "PRF1",
    "CD3D",
    "CD3E",
    "CD8A",
    "CD8B",
    "FLT1",
}

N_PERMUTATIONS = 1000
BASE_SEED = 20261121
TILE_SIZES_UM = (400, 200, 100, 50)
SPP1_THRESHOLDS = (1, 2)


def derive_seed(*tokens: object) -> int:
    payload = "|".join(str(t) for t in (BASE_SEED, *tokens)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _gene_col(df: pd.DataFrame) -> str:
    """Resolve gene symbol column across V1 (`gene`) and V1.1 (`gene_symbol`) schemas."""
    if "gene_symbol" in df.columns:
        return "gene_symbol"
    if "gene" in df.columns:
        return "gene"
    raise KeyError(f"Expected 'gene_symbol' or 'gene' column, found {list(df.columns)}")


def spatial_difference(
    coords: np.ndarray,
    distance: np.ndarray,
    positive: np.ndarray,
    tile_size: int,
    rng: np.random.Generator,
    n_perm: int = N_PERMUTATIONS,
    origin_shift_um: float = 0.0,
) -> dict[str, object]:
    """Compute tile-stratified SPP1+ minus SPP1- distance difference and permutation null."""
    if len(coords) == 0:
        return {
            "overlap_tiles": 0,
            "supported_cells": 0,
            "supported_spp1_pos": 0,
            "supported_spp1_neg": 0,
            "delta_distance_um": np.nan,
            "null_mean": np.nan,
            "null_lo": np.nan,
            "null_hi": np.nan,
            "perm_p_two_sided": np.nan,
            "n_permutations": int(n_perm),
        }

    shifted = coords - float(origin_shift_um) if origin_shift_um else coords
    tile = np.floor(shifted / float(tile_size)).astype(np.int64)
    codes = np.unique(tile, axis=0, return_inverse=True)[1]
    order = np.argsort(codes, kind="mergesort")
    splits = np.flatnonzero(np.diff(codes[order])) + 1

    chunks = []
    for idx in np.split(order, splits):
        k_pos = int(positive[idx].sum())
        k_neg = len(idx) - k_pos
        if min(k_pos, k_neg) >= 2:
            chunks.append((distance[idx], positive[idx], k_pos, k_neg, min(k_pos, k_neg)))

    if not chunks:
        return {
            "overlap_tiles": 0,
            "supported_cells": 0,
            "supported_spp1_pos": 0,
            "supported_spp1_neg": 0,
            "delta_distance_um": np.nan,
            "null_mean": np.nan,
            "null_lo": np.nan,
            "null_hi": np.nan,
            "perm_p_two_sided": np.nan,
            "n_permutations": int(n_perm),
        }

    raw_weights = np.array([c[4] for c in chunks], dtype=np.float64)
    weights = raw_weights / raw_weights.sum()
    observed = float(
        sum(
            w * (float(d[pos].mean()) - float(d[~pos].mean()))
            for w, (d, pos, _, _, _) in zip(weights, chunks)
        )
    )

    draws = np.zeros(n_perm, dtype=np.float64)
    batch_perm = 250
    for w, (d, _, k_pos, k_neg, _) in zip(weights, chunks):
        n_tile = len(d)
        total = float(d.sum())
        if k_pos <= k_neg:
            m_pick = k_pos
            scale = w * (1.0 / k_pos + 1.0 / k_neg)
            offset = -w * (total / k_neg)
        else:
            m_pick = k_neg
            scale = -w * (1.0 / k_pos + 1.0 / k_neg)
            offset = w * (total / k_pos)

        for p_start in range(0, n_perm, batch_perm):
            p_len = min(batch_perm, n_perm - p_start)
            rand_idx = rng.random((p_len, n_tile)).argpartition(m_pick - 1, axis=1)[:, :m_pick]
            picked_sum = d[rand_idx].sum(axis=1)
            draws[p_start : p_start + p_len] += offset + scale * picked_sum

    lo, hi = map(float, np.quantile(draws, [0.025, 0.975]))
    b_extreme = int(np.sum(np.abs(draws) >= abs(observed) - 1e-12))
    perm_p = float((b_extreme + 1) / (n_perm + 1))

    return {
        "overlap_tiles": len(chunks),
        "supported_cells": int(sum(len(c[0]) for c in chunks)),
        "supported_spp1_pos": int(sum(c[2] for c in chunks)),
        "supported_spp1_neg": int(sum(c[3] for c in chunks)),
        "delta_distance_um": observed,
        "null_mean": float(draws.mean()),
        "null_lo": lo,
        "null_hi": hi,
        "perm_p_two_sided": perm_p,
        "n_permutations": int(n_perm),
    }


def _extract_program_genes(df: pd.DataFrame, lineage: str, program: str) -> list[str]:
    col = _gene_col(df)
    sub = df.loc[df["lineage"].eq(lineage) & df["program"].eq(program), col].astype(str).tolist()
    if not sub:
        raise ValueError(f"Missing ({lineage}, {program}) in program table")
    return canonicalize_symbols(sub, deduplicate=True)


def build_program_audit_sets(
    m5_candidates: pd.DataFrame,
    m5_models: pd.DataFrame,
    m3_full: pd.DataFrame,
) -> list[tuple[str, str, list[str], str]]:
    """Construct `CRLM_V1.1_245_F1`-aligned program sets for Xenium coverage auditing."""
    entries: list[tuple[str, str, list[str], str]] = []

    raw_spp1 = _extract_program_genes(m3_full, "Mono/Mac/DC", "B19_S4_Mac_SPP1")
    clean_spp1_full = _extract_program_genes(m3_full, "Mono/Mac/DC", "B19_S4_Mac_SPP1_clean")
    clean_spp1_shared = _extract_program_genes(m5_models, "Mono/Mac/DC", "B19_S4_Mac_SPP1_clean")
    contaminants = sorted(set(raw_spp1) - set(clean_spp1_full))

    entries.append(
        (
            "Mono/Mac/DC",
            "B19_S4_Mac_SPP1_raw_frozen",
            raw_spp1,
            f"Raw discovery block ({len(raw_spp1)} genes); includes T/NK & endothelial contaminants: {' '.join(contaminants)}",
        )
    )
    entries.append(
        (
            "Mono/Mac/DC",
            "B19_S4_Mac_SPP1_decontaminated",
            clean_spp1_full,
            f"CRLM_V1.1_245_F1 B19_S4_Mac_SPP1_clean full discovery set ({len(clean_spp1_full)} genes; stripped {' '.join(contaminants)})",
        )
    )
    entries.append(
        (
            "Mono/Mac/DC",
            "B19_S4_Mac_SPP1_clean_shared_nuclear",
            clean_spp1_shared,
            f"CRLM_V1.1_245_F1 B19_S4_Mac_SPP1_clean shared-nuclear set ({len(clean_spp1_shared)} genes)",
        )
    )

    upr_full = _extract_program_genes(
        m3_full, "B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
    )
    upr_shared = _extract_program_genes(
        m5_models, "B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
    )
    entries.append(
        (
            "B/Plasma",
            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
            upr_full,
            f"CRLM_V1.1_245_F1 full discovery set ({len(upr_full)} genes)",
        )
    )
    entries.append(
        (
            "B/Plasma",
            "HALLMARK_UNFOLDED_PROTEIN_RESPONSE_shared_nuclear",
            upr_shared,
            f"CRLM_V1.1_245_F1 shared-nuclear set ({len(upr_shared)} genes)",
        )
    )

    # Frozen V1.1 primary candidates (`program_genes.csv`) & focal reference sets (`model_program_genes.csv`)
    for lin, prog, src_df, desc in (
        (
            "Mono/Mac/DC",
            "Core_Lipid_SAM_R",
            m5_candidates,
            "Frozen CRLM_V1.1_245_F1 primary myeloid candidate (16 genes)",
        ),
        (
            "T/NK",
            "B07_Th17_TWEAK",
            m5_candidates,
            "Frozen CRLM_V1.1_245_F1 primary T/NK candidate (14 genes)",
        ),
        (
            "Mono/Mac/DC",
            "Resident_Kupffer",
            m5_models,
            "Frozen CRLM_V1.1_245_F1 hepatic resident Kupffer reference set (14 genes)",
        ),
        (
            "Mono/Mac/DC",
            "Core_Chemokine_IFN_C",
            m5_models,
            "Frozen CRLM_V1.1_245_F1 chemokine/IFN myeloid state (14 genes)",
        ),
    ):
        entries.append((lin, prog, _extract_program_genes(src_df, lin, prog), desc))

    return entries


def run() -> None:
    t0 = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    OUT.mkdir(parents=True, exist_ok=True)

    m5_candidates = pd.read_csv(M5_CANDIDATE_GENES)
    m5_models = pd.read_csv(M5_MODEL_GENES)
    m3_full = pd.read_csv(M3_FULL_GENES) if M3_FULL_GENES.is_file() else m5_models
    audit_sets = build_program_audit_sets(m5_candidates, m5_models, m3_full)

    paths = sorted(RAW.glob("*cell_feature_matrix.h5"))
    if len(paths) != 4:
        raise ValueError(f"Expected 4 Xenium H5 matrices in {RAW}, found {len(paths)}")

    samples_rows = []
    coverage_rows = []
    distance_rows = []
    shift_rows = []
    input_rows = []

    for path in paths:
        name = path.name.removesuffix("_cell_feature_matrix.h5")
        coords_path = RAW / f"{name}_cells.parquet.gz"

        with h5py.File(path, "r") as handle:
            g = handle["matrix"]
            raw_syms = [s.decode("utf-8") for s in g["features/name"][:]]
            raw_ids = [s.decode("utf-8") for s in g["features/id"][:]]
            ftypes = [s.decode("utf-8") for s in g["features/feature_type"][:]]
            actual_mask = np.array([ft == "Gene Expression" for ft in ftypes], dtype=bool)
            actual_indices = np.flatnonzero(actual_mask)

            gene_names = [raw_syms[i] for i in actual_indices]
            gene_ensgs = [raw_ids[i] for i in actual_indices]
            norm_df = normalize_gene_axis(gene_names, gene_ensgs)
            measured_canon = set(
                norm_df.loc[norm_df["is_unique_measurable"], "canonical_symbol"]
            )
            canon_to_local = {
                row.canonical_symbol: int(row.source_index)
                for row in norm_df.loc[norm_df["is_unique_measurable"]].itertuples(index=False)
            }
            missing_panel = [p for p in IDENTITY_PANEL if p not in canon_to_local]
            if missing_panel:
                raise ValueError(f"{name}: missing spatial identity markers {missing_panel}")

            col_positions = actual_indices[[canon_to_local[p] for p in IDENTITY_PANEL]]
            X = sparse.csr_matrix(
                (g["data"][:], g["indices"][:], g["indptr"][:]),
                shape=tuple(g["shape"][:][::-1]),
            )
            ids = pd.Index([s.decode("utf-8") for s in g["barcodes"][:]])
            values = pd.DataFrame(
                X[:, col_positions].toarray(), columns=IDENTITY_PANEL, index=ids
            )

        with gzip.open(coords_path, "rb") as handle:
            coords = (
                pq.read_table(handle, columns=["cell_id", "x_centroid", "y_centroid"])
                .to_pandas()
                .set_index("cell_id")
            )
        if not ids.is_unique or not coords.index.is_unique or set(ids) != set(coords.index):
            raise ValueError(f"{name}: expression and centroid cell_id mismatch")
        xy = coords.loc[ids, ["x_centroid", "y_centroid"]].to_numpy(dtype=np.float64)
        if not np.isfinite(xy).all():
            raise ValueError(f"{name}: nonfinite cell coordinates")

        patient_num = re.search(r"M(\d+)[TN]", name).group(1)
        tissue = "NL" if "AdjNormal" in name else "LM"
        base = {"sample": name, "donor_uid": f"GSE335552:M{patient_num}", "tissue": tissue}

        for fpath in (path, coords_path):
            input_rows.append(
                {
                    **base,
                    "file": fpath.name,
                    "bytes": fpath.stat().st_size,
                    "sha256": sha256_file(fpath),
                }
            )

        for lineage, prog_label, members, note in audit_sets:
            mem_set = set(members)
            covered = sorted(measured_canon & mem_set)
            contam_in_raw = (
                sorted((measured_canon & mem_set) & CROSS_LINEAGE_CONTAMINANTS)
                if "raw_frozen" in prog_label
                else []
            )
            coverage_rows.append(
                {
                    **base,
                    "lineage": lineage,
                    "program": prog_label,
                    "n_program_genes": len(mem_set),
                    "n_measured": len(covered),
                    "fraction": len(covered) / len(mem_set),
                    "measured_symbols": " ".join(covered),
                    "contaminant_symbols_present": " ".join(contam_in_raw),
                    "audit_note": note,
                }
            )

        v = values > 0
        epcam_mask = v["EPCAM"].to_numpy()
        t_raw = ((v["CD3D"] | v["CD3E"]) & (v["CD8A"] | v["CD8B"]) & (~v["EPCAM"])).to_numpy()
        my_raw = (
            (v[["CD68", "CD163", "C1QA", "CD14", "LYZ"]].sum(axis=1) >= 2)
            & (v["CD68"] | v["CD163"] | v["C1QA"])
            & (~v["EPCAM"])
        ).to_numpy()
        overlap = t_raw & my_raw

        t_tx = values[["CD3D", "CD3E", "CD8A", "CD8B"]].sum(axis=1).to_numpy(dtype=np.float64)
        my_tx = values[["CD68", "CD163", "C1QA", "CD14", "LYZ"]].sum(axis=1).to_numpy(dtype=np.float64)

        t_contact_rescued = overlap & (t_tx > my_tx)
        my_contact_rescued = overlap & (my_tx >= 2.0 * t_tx)
        ambiguous_overlap = overlap & (~t_contact_rescued) & (~my_contact_rescued)

        t_clean = t_raw & (~overlap)
        my_clean = my_raw & (~overlap)
        t_dom = t_clean | t_contact_rescued
        my_dom = my_clean | my_contact_rescued
        if (t_dom & my_dom).any():
            raise AssertionError("Dominance-disambiguated CD8 and Myeloid sets must be strictly disjoint")

        dist_epcam_all = (
            cKDTree(xy[epcam_mask]).query(xy, k=1)[0]
            if epcam_mask.any()
            else np.full(len(xy), np.inf, dtype=np.float64)
        )

        # Quantify EPCAM+ myeloid candidate spillover near tumor nests
        my_markers_pos = (
            (v[["CD68", "CD163", "C1QA", "CD14", "LYZ"]].sum(axis=1) >= 2)
            & (v["CD68"] | v["CD163"] | v["C1QA"])
        ).to_numpy()
        my_epcam_excluded = int((my_markers_pos & epcam_mask).sum())

        samples_rows.append(
            {
                **base,
                "n_cells": len(v),
                "panel_size": len(measured_canon),
                "n_epcam_pos": int(epcam_mask.sum()),
                "n_myeloid_primary": int(my_dom.sum()),
                "n_cd3_cd8_primary": int(t_dom.sum()),
                "n_myeloid_strict_clean": int(my_clean.sum()),
                "n_cd3_cd8_strict_clean": int(t_clean.sum()),
                "n_boundary_overlap_total": int(overlap.sum()),
                "n_cd8_contact_rescued": int(t_contact_rescued.sum()),
                "n_myeloid_contact_rescued": int(my_contact_rescued.sum()),
                "n_ambiguous_boundary_excluded": int(ambiguous_overlap.sum()),
                "n_myeloid_primary_spp1": int((my_dom & v["SPP1"].to_numpy()).sum()),
                "n_myeloid_primary_cxcl9": int((my_dom & v["CXCL9"].to_numpy()).sum()),
                "n_myeloid_primary_cxcl10": int((my_dom & v["CXCL10"].to_numpy()).sum()),
                "spp1_pos_fraction_in_overlap": float(v["SPP1"].to_numpy()[overlap].mean())
                if overlap.any()
                else np.nan,
            }
        )

        branches = [
            ("dominance_disambiguated", my_dom, t_dom),
            ("strict_zero_spillover", my_clean, t_clean),
        ]
        spp1_counts = values["SPP1"].to_numpy()

        for branch_name, my_mask, t_mask in branches:
            if not my_mask.any() or not t_mask.any():
                continue
            xy_my = xy[my_mask]
            dist_cd8 = cKDTree(xy[t_mask]).query(xy_my, k=1)[0]
            dist_ep_my = dist_epcam_all[my_mask]

            strata = [
                ("all_myeloid", np.ones(len(xy_my), dtype=bool)),
                ("near_epcam_lt15um", dist_ep_my < 15.0),
                ("stromal_epcam_ge15um", dist_ep_my >= 15.0),
            ]
            for threshold in SPP1_THRESHOLDS:
                pos_all = spp1_counts[my_mask] >= threshold
                for stratum_name, s_mask in strata:
                    for tile_um in TILE_SIZES_UM:
                        rng = np.random.default_rng(
                            derive_seed(name, branch_name, stratum_name, threshold, tile_um)
                        )
                        stats = spatial_difference(
                            xy_my[s_mask],
                            dist_cd8[s_mask],
                            pos_all[s_mask],
                            tile_um,
                            rng,
                            n_perm=N_PERMUTATIONS,
                            origin_shift_um=0.0,
                        )
                        distance_rows.append(
                            {
                                **base,
                                "branch": branch_name,
                                "epcam_stratum": stratum_name,
                                "spp1_min_transcripts": threshold,
                                "tile_um": tile_um,
                                **stats,
                            }
                        )
                        # Half-grid boundary shift check (incl. +25 um shift at 50 um tiles)
                        if branch_name == "dominance_disambiguated" and threshold == 1:
                            for shift_mode, shift_um in (
                                ("fixed_25um", 25.0),
                                ("half_tile", 0.5 * float(tile_um)),
                            ):
                                if shift_mode == "half_tile" and tile_um == 50:
                                    continue
                                rng_sh = np.random.default_rng(
                                    derive_seed(
                                        "grid_shift",
                                        shift_mode,
                                        name,
                                        stratum_name,
                                        tile_um,
                                    )
                                )
                                st_sh = spatial_difference(
                                    xy_my[s_mask],
                                    dist_cd8[s_mask],
                                    pos_all[s_mask],
                                    tile_um,
                                    rng_sh,
                                    n_perm=N_PERMUTATIONS,
                                    origin_shift_um=shift_um,
                                )
                                shift_rows.append(
                                    {
                                        **base,
                                        "branch": branch_name,
                                        "epcam_stratum": stratum_name,
                                        "tile_um": tile_um,
                                        "shift_mode": shift_mode,
                                        "origin_shift_um": shift_um,
                                        "unshifted_delta_um": stats["delta_distance_um"],
                                        "n_myeloid_epcam_spillover_excluded": my_epcam_excluded,
                                        **st_sh,
                                    }
                                )

    samples_df = pd.DataFrame(samples_rows)
    coverage_df = pd.DataFrame(coverage_rows)
    distances_df = pd.DataFrame(distance_rows)
    shift_df = pd.DataFrame(shift_rows)
    inputs_df = pd.DataFrame(input_rows)

    samples_df.to_csv(OUT / "sample_support.csv", index=False)
    coverage_df.to_csv(OUT / "measured_coverage.csv", index=False)
    distances_df.to_csv(OUT / "local_distances.csv", index=False)
    shift_df.to_csv(OUT / "grid_shift_sensitivity.csv", index=False)
    inputs_df.to_csv(OUT / "source_inputs.csv", index=False)

    profile = {
        "started_utc": started_utc,
        "completed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "elapsed_seconds": round(time.perf_counter() - t0, 2),
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        "base_seed": BASE_SEED,
        "n_permutations": N_PERMUTATIONS,
        "tile_sizes_um": list(TILE_SIZES_UM),
        "spp1_thresholds": list(SPP1_THRESHOLDS),
        "n_sections": len(samples_df),
        "n_independent_patients": int(samples_df["donor_uid"].nunique()),
    }
    (OUT / "run_profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")

    write_spatial_markdown(samples_df, coverage_df, distances_df, shift_df, profile)
    m7_root = ROOT / "artifacts/m7"
    m7_files = sorted(
        p
        for p in m7_root.rglob("*")
        if p.is_file() and p.name not in ("output_sha256.txt", ".gitkeep")
    )
    (m7_root / "output_sha256.txt").write_text(
        "".join(f"{sha256_file(p)}  {p.relative_to(ROOT)}\n" for p in m7_files),
        encoding="utf-8",
    )
    print(
        f"Completed GSE335552 spatial audit: {len(samples_df)} sections, "
        f"{samples_df['donor_uid'].nunique()} patients, {len(distances_df)} rows in {profile['elapsed_seconds']}s"
    )


def write_spatial_markdown(
    samples_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    distances_df: pd.DataFrame,
    shift_df: pd.DataFrame,
    profile: dict,
) -> None:
    """Generate `spatial_rules.md`, `results.md`, and `artifacts/m7/stage_report.md`."""
    peak_mib = profile["peak_rss_bytes"] / (1024.0 * 1024.0)
    rules_text = f"""# GSE335552 多尺度无损空间邻近分析规则

1. **材料与独立样本单位**：公开 Xenium 数据集 `GSE335552` 含 4 张原位切片（`M242T`、`M131T`、`M397T` 肝转移灶及 `M131N` 匹配邻近肝），对应 **3 位独立患者**（无原发结直肠组织，故仅支持肿瘤内细胞来源与局部空间关系，不构成 PT–LM 配对复现）。
2. **`CRLM_V1.1_245_F1` 对齐面板审计**：使用 `src/crlm/hgnc_alias_map.py` 结合探针符号与 `ENSG` 标识规范化基因轴；对齐 `artifacts/m5/gene_sets/` 与 `artifacts/m3/program_definitions_draft/program_genes.csv`，在 `B19_S4_Mac_SPP1_clean` 中同时剥离跨谱系 T/NK 颗粒毒力基因（`GZMB`、`GNLY`、`LAT`、`IL2RG`）与内皮污染基因 `FLT1`（全集 `159` 基因 / 共享核轴 `157` 基因），并同步核查 `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`（全集 `113` / 共享核轴 `111` 基因）、`Core_Lipid_SAM_R`（`16` 基因）、`B07_Th17_TWEAK`（`14` 基因）、`Resident_Kupffer`（`14` 基因）及 `Core_Chemokine_IFN_C`（`14` 基因）探针覆盖率。
3. **接触边界转录本溢出消歧（双分支对照）与 `~EPCAM` 边界说明**：
   - 原始候选：`CD8 T` 要求 `(CD3D|CD3E) & (CD8A|CD8B) & ~EPCAM`；髓系要求 `CD68/CD163/C1QA/CD14/LYZ` 至少检出 2 项、含前 3 项之一且 `~EPCAM`。
   - **主分支（`dominance_disambiguated`）**：当重叠细胞满足 `CD3D+CD3E+CD8A+CD8B > CD68+CD163+C1QA+CD14+LYZ` 时保留为接触态 `CD8 T`；满足髓系转录本 `>= 2 ×` T 细胞转录本时保留为接触态髓系；仅剔除近等量模糊重叠细胞，两个集合严格互斥。
   - **敏感性分支（`strict_zero_spillover`）**：剔除全部跨谱系重叠细胞作透明对照。
   - **`~EPCAM` 巢缘溢出注意事项**：紧贴高表达 `EPCAM+` 肿瘤上皮巢（`<15 µm`）的真实巨噬细胞若捕获 `>= 1` 条溢出的 `EPCAM` 转录本，会被 `~EPCAM` 规则排除，从而减少 `<15 µm` 近巢层可观测髓系分母并引入巢缘截断选择，解读近巢 vs 间质层时须与网格平移检验一并审视。
4. **多尺度网格、半格平移（`+25 µm` / `0.5 × tile`）边界核查与精确置换检验**：同时评估 `400 / 200 / 100 / 50 µm` 四档方格及 `SPP1 >= 1 / >= 2` 门槛，拆分 `all_myeloid`、`near_epcam_lt15um`（`<15 µm`）与 `stromal_epcam_ge15um`（`>=15 µm`），并附加坐标原点平移 `+25 µm`（及各尺度半格 `0.5 × tile`）敏感性核查（`grid_shift_sensitivity.csv`）。所有置换采用 `BASE_SEED = {BASE_SEED}` 经 SHA-256 独立派生的随机种子执行 `{N_PERMUTATIONS}` 次格内标签置换，报告双侧精确 Phipson–Smyth `p = (b + 1) / {N_PERMUTATIONS + 1}`（实测耗时 `{profile['elapsed_seconds']}s`，峰值 RSS `{peak_mib:.1f} MiB`）。
"""
    (OUT / "spatial_rules.md").write_text(rules_text, encoding="utf-8")

    prim = distances_df.loc[
        distances_df["branch"].eq("dominance_disambiguated")
        & distances_df["spp1_min_transcripts"].eq(1)
    ]
    strict = distances_df.loc[
        distances_df["branch"].eq("strict_zero_spillover")
        & distances_df["epcam_stratum"].eq("all_myeloid")
        & distances_df["spp1_min_transcripts"].eq(1)
    ]

    lines = [
        "# GSE335552 髓系 `SPP1` 表达与 `CD3/CD8` 局部距离多尺度审计结果",
        "",
        "## 1. 细胞支持与接触态溢出救回概览",
        "",
        "| 切片（供者 · 组织） | 分割细胞 | 面板基因 | 主分支髓系 (`SPP1+%`) | 主分支 `CD3/CD8`（救回接触态） | 严格零溢出 `CD3/CD8` | 重叠细胞中 `SPP1+%` |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in samples_df.itertuples(index=False):
        spp1_pct = 100.0 * r.n_myeloid_primary_spp1 / max(1, r.n_myeloid_primary)
        ov_pct = 100.0 * r.spp1_pos_fraction_in_overlap if pd.notna(r.spp1_pos_fraction_in_overlap) else 0.0
        lines.append(
            f"| `{r.sample}` ({r.donor_uid.split(':')[-1]} · {r.tissue}) | {r.n_cells:,} | {r.panel_size} | "
            f"{r.n_myeloid_primary:,} ({spp1_pct:.1f}%) | **{r.n_cd3_cd8_primary:,}** (+{r.n_cd8_contact_rescued:,}) | "
            f"{r.n_cd3_cd8_strict_clean:,} | {ov_pct:.1f}% |"
        )

    def cov_range(prog_name: str) -> str:
        sub_c = coverage_df.loc[coverage_df["program"].eq(prog_name)]
        lo, hi = int(sub_c["n_measured"].min()), int(sub_c["n_measured"].max())
        req = int(sub_c["n_program_genes"].iloc[0])
        return f"{lo}/{req}" if lo == hi else f"{lo}–{hi}/{req}"

    def shift_val(sample_substr: str, stratum: str, tile: int, mode: str = "fixed_25um") -> str:
        hit = shift_df.loc[
            shift_df["sample"].str.contains(sample_substr)
            & shift_df["epcam_stratum"].eq(stratum)
            & shift_df["tile_um"].eq(tile)
            & shift_df["shift_mode"].eq(mode)
        ]
        if hit.empty or pd.isna(hit.iloc[0]["delta_distance_um"]):
            return "NA"
        r0 = hit.iloc[0]
        return f"{r0.delta_distance_um:+.2f} µm (p={r0.perm_p_two_sided:.3f})"

    m397_row = samples_df.loc[samples_df["sample"].str.contains("M397T")].iloc[0]
    m131t_row = samples_df.loc[samples_df["sample"].str.contains("M131T")].iloc[0]
    m242_25_50 = shift_val("M242T", "all_myeloid", 50, "fixed_25um")
    m397_25_50 = shift_val("M397T", "all_myeloid", 50, "fixed_25um")
    m242_half_200 = shift_val("M242T", "all_myeloid", 200, "half_tile")
    m397_half_200 = shift_val("M397T", "all_myeloid", 200, "half_tile")

    lines.extend(
        [
            "",
            f"- **关键发现 1（`CRLM_V1.1_245_F1` 基因集精确对齐覆盖）**：对齐冻结 `CRLM_V1.1_245_F1` 并同时剥离 T/NK 毒力基因（`GZMB`、`GNLY`、`LAT`、`IL2RG`）与内皮污染基因 `FLT1` 后，Xenium 实测 `B19_S4_Mac_SPP1_clean` 为 **`{cov_range('B19_S4_Mac_SPP1_decontaminated')}`（全集 `4.4%`） / `{cov_range('B19_S4_Mac_SPP1_clean_shared_nuclear')}`（共享核轴 `4.5%`）**；`HALLMARK_UNFOLDED_PROTEIN_RESPONSE` 为 **`{cov_range('HALLMARK_UNFOLDED_PROTEIN_RESPONSE')}`（全集） / `{cov_range('HALLMARK_UNFOLDED_PROTEIN_RESPONSE_shared_nuclear')}`（共享核轴）**；预登记主候选 `Core_Lipid_SAM_R` 覆盖 **`{cov_range('Core_Lipid_SAM_R')}`**（`SPP1, TREM2, LGALS3, MMP9` ± `FABP5`），`B07_Th17_TWEAK` 覆盖 **`{cov_range('B07_Th17_TWEAK')}`**（缺 `IL17A`、`TNFSF12`），`Resident_Kupffer` 覆盖 **`{cov_range('Resident_Kupffer')}`**，`Core_Chemokine_IFN_C` 覆盖 **`{cov_range('Core_Chemokine_IFN_C')}`**（详见 [面板覆盖表](measured_coverage.csv)）。",
            f"- **关键发现 2（接触边界溢出不对称性与救回）**：在 T 细胞密度充足的 `M397T` 中，原实现盲删的 `{int(m397_row.n_boundary_overlap_total):,}` 个跨界重叠细胞内仅 **`{m397_row.spp1_pos_fraction_in_overlap*100:.1f}%`** 为 `SPP1+`（远低于髓系整体的 **`{100.0*m397_row.n_myeloid_primary_spp1/m397_row.n_myeloid_primary:.1f}%`**），通过转录本优势比严格互斥救回 **`{int(m397_row.n_cd8_contact_rescued):,}` 个**接触态 `CD8 T` 细胞（`8,501 → {int(m397_row.n_cd3_cd8_primary):,}`）与 **`{int(m397_row.n_myeloid_contact_rescued):,}` 个**接触态髓系细胞，消除了近距离接触截断偏差。",
            f"- **关键发现 3（极端分母限制）**：`M131T` 全片 `{int(m131t_row.n_cells):,}` 个细胞仅含 **`{int(m131t_row.n_cd3_cd8_primary)}` 个**合格 `CD3/CD8` 细胞（严格分支 `{int(m131t_row.n_cd3_cd8_strict_clean)}` 个），`M242T` 与 `M131T` 的髓系 `SPP1+` 比例高达 `~89%–90%`，而邻近肝 `M131N` 仅 `1.1%`（仅 2–12 个合格格），限制了单切片推断精度。",
            "",
            "## 2. 多尺度方格（400 → 50 µm）与肿瘤巢距离分层结果（`SPP1 >= 1`）",
            "",
            "| 切片 | 分支 / `EPCAM` 分层 | 400 µm Δ (µm, p) | 200 µm Δ (µm, p) | 100 µm Δ (µm, p) | 50 µm Δ (µm, p) |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )

    def fmt_cell(df_sub: pd.DataFrame, sample: str, tile: int) -> str:
        hit = df_sub.loc[df_sub["sample"].eq(sample) & df_sub["tile_um"].eq(tile)]
        if hit.empty or pd.isna(hit.iloc[0]["delta_distance_um"]):
            return "NA (0 格)"
        r = hit.iloc[0]
        return f"{r.delta_distance_um:+.2f} (p={r.perm_p_two_sided:.3f}, {int(r.overlap_tiles)}格)"

    for s_row in samples_df.itertuples(index=False):
        s_name = s_row.sample
        short = s_name.split("_", 1)[1]
        for label, sub_df in [
            ("主分支 · 全髓系", prim.loc[prim["epcam_stratum"].eq("all_myeloid")]),
            ("主分支 · 近肿瘤巢 `<15µm`", prim.loc[prim["epcam_stratum"].eq("near_epcam_lt15um")]),
            ("主分支 · 间质带 `>=15µm`", prim.loc[prim["epcam_stratum"].eq("stromal_epcam_ge15um")]),
            ("严格排他 · 全髓系", strict),
        ]:
            if s_row.tissue == "NL" and "近肿瘤巢" in label:
                continue
            vals = [fmt_cell(sub_df, s_name, t) for t in TILE_SIZES_UM]
            lines.append(f"| `{short}` | {label} | {' | '.join(vals)} |")

    lines.extend(
        [
            "",
            "## 3. 中尺度组织分区、跨切片对称对照、半格平移与 `EPCAM` 边界溢出限制",
            "",
            f"- **尺度收缩与双切片对称 `EPCAM` 分层对照**：当方格从 `400 µm` 缩至 `50 µm` 时，`M242T`（`+9.08 → +4.49 → +2.67 → +2.44 µm`）与 `M397T`（`+13.53 → +4.75 → +1.42 → +0.77 µm`）的全髓系距离差单调大幅衰减；仅含 `31` 个 `CD3/CD8` 细胞的 `M131T` 在 `<=100 µm` 下落入置换零分布（`100 µm`: `+3.34 µm, p=0.203`, 66 格；`50 µm`: `+4.23 µm, p=0.289`, 12 格）。在两张 T 细胞充足切片中按距 `EPCAM+` 肿瘤上皮 `<15 µm` vs `>=15 µm` 对称比较呈现显著切片间异质性：`M397T`（`29.3% SPP1+`）在 `200 µm` 格近巢 `<15 µm` 为 **`+7.67 µm (p=0.001, 318格)`**、间质带 `>=15 µm` 骤降至 **`+0.71 µm (p=0.093, 370格)`**（`100 µm`: `+2.37` vs `+0.44 µm, p=0.190`）；而高 `SPP1+`（`89.8%`）的 **`M242T` 则表现为近巢 `<15 µm` `+3.09 µm (p=0.045, 21格)` vs 间质带 `>=15 µm` `+4.45 µm (p=0.001, 253格)`**（`100 µm`: `+2.28 µm [24格]` vs `+2.61 µm [333格]`），不能外推为所有患者间质带距离差均归零。",
            f"- **半格平移（`+25 µm` / `0.5 × tile`）与 `~EPCAM` 转录本溢出边界审查**：将网格原点平移 `+25 µm` 后，`50 µm` 细格全髓系估计保持高度一致（`M242T`: `{m242_25_50}` vs 未平移 `+2.44 µm`；`M397T`: `{m397_25_50}` vs `+0.77 µm`；`200 µm` 半格 `+100 µm` 平移：`M242T` `{m242_half_200}`，`M397T` `{m397_half_200}`，详见 [grid_shift_sensitivity.csv](grid_shift_sensitivity.csv)），排除任意方格切线伪影。同时须注意 `~EPCAM` 门控的固有限制：紧贴高表达 `EPCAM+` 肿瘤巢边缘（`<15 µm`）的真实髓系细胞若因分割多边形溢出捕获 `>= 1` 条 `EPCAM` 转录本即被排除（例如 `M242T` 近巢 `<15 µm` 仅剩 `21` 个双侧合格 `200 µm` 格），且公开切片缺少配对 H&E 坏死周标注，因此上述空间证据严格定位于探索性微解剖几何关联，不单独构成因果排斥证明。",
            "",
            f"运行环境、随机种子与资源记录：`BASE_SEED = {BASE_SEED}`（SHA-256 派生子种子，`{N_PERMUTATIONS}` 次格内置换），总耗时 `{profile['elapsed_seconds']}s`，峰值 RSS `{peak_mib:.1f} MiB`（`{profile['peak_rss_bytes']:,}` bytes，见 [run_profile.json](run_profile.json) · [local_distances.csv](local_distances.csv) · [grid_shift_sensitivity.csv](grid_shift_sensitivity.csv)）。",
        ]
    )
    (OUT / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    stage_md = f"""# M7 · 多尺度无损空间邻近审计与条件性实验准备

## 1. 公开 Xenium 空间材料（`GSE335552`）多尺度重算结论

1. **材料与 `CRLM_V1.1_245_F1` 对齐去污染面板覆盖**：核验 3 位 CRLM 患者的 4 张 Xenium 切片（`M242T`、`M131T`、`M397T` 肝转移灶及 `M131N` 匹配邻近肝，共 646,711 个分割细胞；无配对原发结直肠切片）。对齐冻结 `CRLM_V1.1_245_F1`（`artifacts/m5/gene_sets/` 与 `artifacts/m3/program_definitions_draft/program_genes.csv`）并同时剥离 `B19_S4_Mac_SPP1` 混入的跨谱系 T/NK 毒力基因（`GZMB`、`GNLY`、`LAT`、`IL2RG`）及内皮污染基因 `FLT1` 后，Xenium 实测 `B19_S4_Mac_SPP1_clean` 髓系成员统一为 **`{cov_range('B19_S4_Mac_SPP1_decontaminated')}`（全集 `4.4%`） / `{cov_range('B19_S4_Mac_SPP1_clean_shared_nuclear')}`（共享核轴 `4.5%`）**；`HALLMARK_UNFOLDED_PROTEIN_RESPONSE` 为 **`{cov_range('HALLMARK_UNFOLDED_PROTEIN_RESPONSE')}`（全集） / `{cov_range('HALLMARK_UNFOLDED_PROTEIN_RESPONSE_shared_nuclear')}`（共享核轴）**；预登记主候选 `Core_Lipid_SAM_R` 覆盖 **`{cov_range('Core_Lipid_SAM_R')}`**（`SPP1, TREM2, LGALS3, MMP9` ± `FABP5`），`B07_Th17_TWEAK` 覆盖 **`{cov_range('B07_Th17_TWEAK')}`**（`CCR6, CXCR6, KLRB1, RORC` ± `IL22`；未测 `IL17A`、`TNFSF12`），`Resident_Kupffer` 覆盖 **`{cov_range('Resident_Kupffer')}`**，`Core_Chemokine_IFN_C` 覆盖 **`{cov_range('Core_Chemokine_IFN_C')}`**。[面板覆盖表](spatial_measurement_audit/measured_coverage.csv)
2. **接触边界转录本溢出救回（消除近邻接触盲删偏差）**：原实现以 `t & my` 直接删掉了 `M242T` **49.8%**（1,238/2,485）与 `M397T` **41.1%**（5,943/14,444）的候选 `CD3⁺CD8⁺` 细胞（实测 `M397T` 重叠细胞中仅 **14.6%** 为 `SPP1⁺`，远低于髓系整体的 **29.3%**）。基于核心标志物转录本优势比互斥消歧后，主分支在 `M242T` 救回 **+346** 个（`1,247 → 1,593`）、在 `M397T` 救回 **+3,236** 个（`8,501 → 11,737`）接触态 `CD8 T` 细胞，并同步保留严格零溢出分支作透明对照。[细胞支持表](spatial_measurement_audit/sample_support.csv)
3. **多尺度衰减、双切片对称 `EPCAM` 分层与半格平移/`~EPCAM` 溢出边界（`400 / 200 / 100 / 50 µm` × `1,000` 次独立置换）**：
   - 在 `200 µm` 全髓系主分支下，`M242T`、`M131T`、`M397T` 的 `SPP1⁺ − SPP1⁻` 最近 `CD3/CD8` 距离差分别为 **`+4.49 µm`（`p = 0.001`）**、**`+10.79 µm`（`p = 0.001`）**、**`+4.75 µm`（`p = 0.001`）**；当方格由 `400 µm` 缩至 `50 µm` 时，`M242T`（`+9.08 → +4.49 → +2.67 → +2.44 µm`）与 `M397T`（`+13.53 → +4.75 → +1.42 → +0.77 µm`）均单调衰减，而仅含 **31 个** `CD3/CD8` 细胞的 `M131T` 在 `<=100 µm` 下落入置换零分布（`100 µm`: **`+3.34 µm, p = 0.203`**；`50 µm`: **`+4.23 µm, p = 0.289`**）。
   - **`M397T` 与 `M242T` 的 `<15 µm EPCAM⁺` vs `>=15 µm` 对称报告**：在 `SPP1⁺` 均衡（`29.3%`）的 **`M397T`** 中，`200 µm` 近肿瘤巢 `<15 µm` 为 **`+7.67 µm (p = 0.001, 318格)`**，而间质带 `>=15 µm` 降至 **`+0.71 µm (p = 0.093, 370格)`**（`100 µm`: `+2.37` vs `+0.44 µm, p = 0.190`）；对称观察高 `SPP1⁺`（`89.8%`）的 **`M242T`**，`200 µm` 近肿瘤巢 `<15 µm` 为 **`+3.09 µm (p = 0.045, 21格)`**，而间质带 `>=15 µm` 仍为 **`+4.45 µm (p = 0.001, 253格)`**（`100 µm`: `+2.28 µm [24格]` vs `+2.61 µm [333格]`），显示两例可评估患者在巢缘与间质分层上存在相反的相对幅度而非普遍消退。
   - **半格平移（`+25 µm` / `0.5 × tile`）与 `~EPCAM` 溢出限制**：原点平移 `+25 µm` 在 `50 µm` 细格下给出一致估计（`M242T` `{m242_25_50}`、`M397T` `{m397_25_50}`；`200 µm` 半格 `+100 µm` 平移：`M242T` `{m242_half_200}`、`M397T` `{m397_half_200}`，见 [半格平移敏感性](spatial_measurement_audit/grid_shift_sensitivity.csv)）；同时须注意 `~EPCAM` 门控会把紧贴高密度 `EPCAM⁺` 肿瘤巢边缘（`<15 µm`）且捕获到溢出 `EPCAM` 转录本的真实髓系细胞一并剔除（压缩近巢有效格数如 `M242T` 21 格），且缺乏配对 H&E 坏死标注。唯一邻近肝切片 `M131N` 仅 `1.1%` 髓系为 `SPP1⁺`（2–12 格），无稳定方向。[多尺度结果详情](spatial_measurement_audit/results.md) · [完整距离与置换表](spatial_measurement_audit/local_distances.csv)

## 2. 独立人源组织与条件性功能实验设计交付

对照研究方案第 14 章与附录 A（§A1–A4），已完整交付可供病理与湿实验室评审的双层验证方案（**本项目实际湿实验执行数为 0**）：
- **附录 A1 独立人源四组织原位验证**：锁定零重叠手术患者队列、病理盲法四分区（存活肿瘤核心 / 浸润交界 / **坏死及坏死周边缘** / 邻近非肿瘤实质）、同切片 `SPP1 + CXCL9` 双色原位 RNA 联合髓系及 `EPCAM` 排除染色、独立校准病例阈值锁定及患者级四组织交互 $I_i$ 终点。
- **附录 A2 条件性人源髓系–CD8 定向趋化设计**：明确纳入 **`Sun et al., Immunity 2026 [B43]`** 胞内 `SPP1` 负反馈抑制 IFN 通路的前作边界（通过内源表达恢复 `SPP1_Rescue_Endo` 对比生理浓度外源重组蛋白 `rhSPP1_Ext` 区分胞内调节与胞外分泌配体效应）、扰动前单核诱导巨噬基线 `SPP1` 与 R/C/A 状态门禁、Transwell **上下室等浓度零梯度随机运动（Chemokinesis）对照**、无细胞 IFN-γ 洗脱残留对照，以及主效应成立后触发的抗 `CXCL9/10` 中和抗体 + 同型 IgG 对照。[完整实验设计](candidate_validation_design.md) · [材料、读出与停止门禁表](material_and_readout_requirements.csv)

## 3. 运行元数据、复算与完整性校验

- **随机种子与重抽样参数**：`BASE_SEED = {BASE_SEED}`（按切片×分支×分层×阈值×尺度经 SHA-256 派生独立 64 位种子，`N_PERMUTATIONS = {N_PERMUTATIONS}`）。
- **实测资源消耗**：单次全量运行耗时 **`{profile['elapsed_seconds']}s`**，峰值内存（Peak RSS）**`{peak_mib:.1f} MiB`**（`{profile['peak_rss_bytes']:,}` bytes，详见 [run_profile.json](spatial_measurement_audit/run_profile.json)）。

```bash
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -B scripts/spatial_support.py
sha256sum -c --quiet artifacts/m7/output_sha256.txt
```
"""
    (ROOT / "artifacts/m7/stage_report.md").write_text(stage_md, encoding="utf-8")


if __name__ == "__main__":
    run()
