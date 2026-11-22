"""Build V1.1 patient program, classifier benchmark, multi-scale spatial, and external transfer figures.

Supports:
- Full production execution (`python scripts/build_paper_figures.py`), writing publication PNG/PDFs
  to `artifacts/m8/figures/` and exact plotted tables to `artifacts/m8/source_data/`.
- Pure in-memory verification (`python scripts/build_paper_figures.py --self-check`), rendering all
  four figure families onto `io.BytesIO` buffers from schema-exact synthetic tables with zero disk writes.
"""

from __future__ import annotations

import argparse
import hashlib
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "artifacts"
OUT = A / "m8"
FIG = OUT / "figures"
DATA = OUT / "source_data"

ARCHIVE_V1_EXT_METRICS = (
    A / "archive/CRLM_V1_245_F1/artifacts/m6/GSE164522/model_predictions/patient_metrics.csv"
)
TECHFIX_EXT_METRICS = A / "techfix_v1_f1/models/techfix_external_classifier_metrics.csv"

PALETTE = ["#175b75", "#b14925", "#577644", "#6b4c9a"]
PLANNED_COHORTS = (
    ("GSE245552", "245\n(Dev)"),
    ("GSE164522", "164\n(Pri)"),
    ("GSE203608", "203608\n(Sec)"),
    ("GSE178318", "178\n(Zero-Exp)"),
)

SINGLE_LINEAGE_METHODS = (
    ("program_elastic", "Program Elastic"),
    ("composition_elastic", "Composition Elastic"),
    ("QC_elastic", "QC Elastic (ctrl)"),
    ("PCA_mean", "PCA mean"),
    ("PCA_subtype_balanced", "PCA subtype-balanced"),
    ("PCA_remove_organ_markers", "PCA − organ markers"),
    ("PCA_attention", "PCA attention"),
    ("PCA32_CTA", "PCA CTA"),
    ("PCA32_HA", "PCA HA"),
    ("PCA32_MixMIL", "PCA MixMIL"),
    ("scGPT_mean", "scGPT mean (512D)"),
    ("scGPT_attention", "scGPT attention (512D)"),
    ("scGPT32_mean", "scGPT32 mean"),
    ("scGPT32_CTA", "scGPT32 CTA"),
    ("scGPT32_HA", "scGPT32 HA"),
    ("scGPT32_MixMIL", "scGPT32 MixMIL"),
    ("scGPT512_CTA", "scGPT512 CTA"),
    ("scGPT512_HA", "scGPT512 HA"),
)

JOINT_LINEAGE_METHODS = (
    ("PCA_mean", "PCA mean"),
    ("PCA_attention", "PCA attention"),
    ("PCA32_CTA", "PCA CTA"),
    ("PCA32_HA", "PCA HA"),
    ("PCA32_MixMIL", "PCA MixMIL"),
    ("scGPT32_mean", "scGPT32 mean"),
    ("scGPT32_attention", "scGPT32 attention"),
    ("scGPT32_CTA", "scGPT32 CTA"),
    ("scGPT32_HA", "scGPT32 HA"),
    ("scGPT32_MixMIL", "scGPT32 MixMIL"),
    ("scGPT_mean", "scGPT mean (512D)"),
    ("scGPT_attention", "scGPT attention (512D)"),
    ("scGPT512_CTA", "scGPT512 CTA"),
    ("scGPT512_HA", "scGPT512 HA"),
)

EXTERNAL_METHODS = (
    ("program_elastic", "Program"),
    ("PCA_mean", "PCA mean"),
    ("PCA_attention", "PCA attention"),
    ("scGPT_mean", "scGPT mean"),
    ("scGPT_attention", "scGPT attention"),
)


def apply_rcparams() -> None:
    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _savefig_metadata(ext: str) -> dict[str, object]:
    if ext == "pdf":
        return {"CreationDate": None, "ModDate": None, "Producer": "Matplotlib"}
    return {"Software": "Matplotlib"}


def save_or_verify(fig: plt.Figure, stem: str, fig_dir: Path | None) -> dict[str, int]:
    """Save PNG + PDF to `fig_dir`, or verify non-empty in-memory rendering when `fig_dir is None`."""
    sizes: dict[str, int] = {}
    try:
        if fig_dir is not None:
            for ext in ("png", "pdf"):
                out_path = fig_dir / f"{stem}.{ext}"
                fig.savefig(
                    out_path,
                    format=ext,
                    dpi=200,
                    bbox_inches="tight",
                    metadata=_savefig_metadata(ext),
                )
                sizes[ext] = out_path.stat().st_size
        else:
            for ext in ("png", "pdf"):
                buf = io.BytesIO()
                fig.savefig(
                    buf,
                    format=ext,
                    dpi=160,
                    bbox_inches="tight",
                    metadata=_savefig_metadata(ext),
                )
                data = buf.getvalue()
                if len(data) < 1500:
                    raise RuntimeError(f"In-memory {stem}.{ext} buffer unexpectedly small ({len(data)} B)")
                sizes[ext] = len(data)
    finally:
        plt.close(fig)
    return sizes


def _load_primary_candidates(candidates_df: pd.DataFrame | None) -> list[tuple[str, str, str]]:
    """Return `[(lineage, program, expected_direction), ...]` for frozen primary candidates."""
    if candidates_df is not None and not candidates_df.empty:
        df = candidates_df.copy()
        if "candidate_tier" in df.columns:
            df = df.loc[df["candidate_tier"].astype(str).eq("primary")]
        if not df.empty:
            out = []
            for r in df.itertuples(index=False):
                exp = str(getattr(r, "expected_direction_DT", getattr(r, "expected_direction", "positive")))
                out.append((str(r.lineage), str(r.program), exp))
            return out
    return [
        ("Mono/Mac/DC", "Core_Lipid_SAM_R", "positive"),
        ("T/NK", "B07_Th17_TWEAK", "negative"),
    ]


def _short_lineage(lineage: str) -> str:
    return {"Mono/Mac/DC": "Myeloid"}.get(lineage, lineage)


def _scatter_with_mean(
    ax: plt.Axes,
    x: float,
    values: np.ndarray,
    color: str,
    *,
    annotate_signs: bool = True,
) -> None:
    n = len(values)
    if n == 0:
        return
    offset = np.linspace(-0.13, 0.13, n) if n > 1 else np.zeros(1)
    ax.scatter(
        x + offset,
        values,
        s=25,
        alpha=0.75,
        color=color,
        edgecolors="none",
        zorder=3,
    )
    mean_v = float(np.mean(values))
    ax.plot([x - 0.22, x + 0.22], [mean_v, mean_v], color="#111", lw=2.0, zorder=4)
    if annotate_signs:
        n_pos = int(np.sum(values > 0))
        n_neg = int(np.sum(values < 0))
        va = "bottom" if mean_v >= 0 else "top"
        y_pos = float(np.max(values)) if mean_v >= 0 else float(np.min(values))
        span = max(abs(float(np.max(values)) - float(np.min(values))), 0.02)
        dy = 0.06 * span if mean_v >= 0 else -0.06 * span
        ax.text(
            x,
            y_pos + dy,
            f"{n_pos}+/{n_neg}−",
            ha="center",
            va=va,
            fontsize=7.2,
            color="#333",
        )


def render_patient_programs(
    patients_df: pd.DataFrame,
    candidates_df: pd.DataFrame | None,
    m3_decomp_df: pd.DataFrame | None,
    m6_decomp_df: pd.DataFrame | None,
    m6_178_decomp_df: pd.DataFrame | None = None,
    *,
    fig_dir: Path | None = FIG,
    data_dir: Path | None = DATA,
) -> dict[str, int]:
    """Build `01_patient_programs` (`2 x 3`: 2 V1.1 primary candidates + UPR falsification control)."""
    primary = _load_primary_candidates(candidates_df)
    target_keys = set((lin, prog) for lin, prog, _ in primary) | {
        ("B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
    }
    keep_mask = [
        (str(lin), str(prog)) in target_keys
        for lin, prog in zip(patients_df["lineage"], patients_df["program"])
    ]
    sub_patients = patients_df.loc[keep_mask].copy()
    if data_dir is not None:
        sub_patients.to_csv(data_dir / "candidate_patient_effects.csv", index=False)

    n_cols = len(primary) + 1
    fig, axes = plt.subplots(2, n_cols, figsize=(4.55 * n_cols, 7.1), layout="constrained")

    # Columns 0..k-1: V1.1 Primary Candidates across all 4 planned cohorts (incl. explicit GSE178318 slot) + Quartet Q
    for col, (lineage, program, exp_dir) in enumerate(primary):
        color = PALETTE[col % len(PALETTE)]
        dir_tag = "LM > PT" if exp_dir.startswith("pos") or exp_dir == "+1" else "LM < PT"
        title_base = f"{_short_lineage(lineage)} · {program}"
        pp = sub_patients.loc[
            sub_patients["lineage"].eq(lineage) & sub_patients["program"].eq(program)
        ]

        # Row 0: Paired D_T across 245, 164, 203608, and explicit 178 placeholder
        ax_top = axes[0, col]
        ax_top.axhline(0, color="#9ca3af", lw=0.85, ls="--", zorder=1)
        xtick_labels = []
        observed_all = pp["D_T"].dropna().to_numpy(dtype=float)
        ymax_abs = max(float(np.max(np.abs(observed_all))) if len(observed_all) else 0.08, 0.04)

        for x, (acc, short_name) in enumerate(PLANNED_COHORTS):
            vals = pp.loc[pp["dataset"].eq(acc), "D_T"].dropna().to_numpy(dtype=float)
            if acc == "GSE178318" and len(vals) == 0:
                ax_top.axvspan(x - 0.34, x + 0.34, color="#f3f4f6", alpha=0.95, zorder=0)
                ax_top.axvline(x - 0.34, color="#cbd5e1", lw=0.8, ls=":")
                ax_top.axvline(x + 0.34, color="#cbd5e1", lw=0.8, ls=":")
                ax_top.text(
                    x,
                    0.0,
                    "Terms pending\n(N=6 unexposed;\nraw p=NA, p_H=1)",
                    ha="center",
                    va="center",
                    fontsize=6.8,
                    color="#4b5563",
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cbd5e1", lw=0.7),
                    zorder=5,
                )
                xtick_labels.append(f"{short_name}\nN = NA")
            else:
                if acc == "GSE178318":
                    ax_top.axvspan(x - 0.35, x + 0.35, color="#ecfdf5", alpha=0.88, zorder=0)
                    ax_top.axvline(x - 0.35, color="#6ee7b7", lw=0.8, ls=":")
                    ax_top.axvline(x + 0.35, color="#6ee7b7", lw=0.8, ls=":")
                _scatter_with_mean(ax_top, x, vals, color)
                xtick_labels.append(f"{short_name}\nN = {len(vals)}")

        ax_top.set_xlim(-0.52, len(PLANNED_COHORTS) - 0.48)
        ax_top.set_ylim(-1.25 * ymax_abs, 1.28 * ymax_abs)
        ax_top.set_xticks(range(len(PLANNED_COHORTS)), xtick_labels, fontsize=7.8)
        ax_top.set_ylabel("Debiased UCell difference (LM − PT)")
        ax_top.set_title(f"{title_base} (exp: {dir_tag})", loc="left", fontsize=9.8, weight="bold")
        ax_top.set_xlabel("Patient-paired D_T across planned cohorts")

        # Row 1: GSE164522 four-tissue quartet Q (D_T_Q, D_N_Q, I)
        ax_bot = axes[1, col]
        ax_bot.axhline(0, color="#9ca3af", lw=0.85, ls="--", zorder=1)
        q = pp.loc[pp["dataset"].eq("GSE164522") & pp["I"].notna()].sort_values("donor_uid")
        if not q.empty:
            for r in q.itertuples(index=False):
                ax_bot.plot(
                    [0, 1, 2],
                    [float(r.D_T_Q), float(r.D_N_Q), float(r.I)],
                    color=color,
                    alpha=0.18,
                    lw=0.85,
                    zorder=2,
                )
        q_labels = []
        for x, (c_col, c_lbl) in enumerate(
            (("D_T_Q", "LM − PT\n(D_T_Q)"), ("D_N_Q", "NL − NC\n(D_N_Q)"), ("I", "Interaction\n(I = D_T − D_N)"))
        ):
            vals = q[c_col].dropna().to_numpy(dtype=float) if not q.empty else np.array([], dtype=float)
            _scatter_with_mean(ax_bot, x, vals, color)
            q_labels.append(f"{c_lbl}\nN = {len(vals)}")
        q_vals_all = (
            q[["D_T_Q", "D_N_Q", "I"]].to_numpy(dtype=float).ravel()
            if not q.empty
            else np.array([0.05])
        )
        q_max = max(float(np.nanmax(np.abs(q_vals_all))), 0.04)
        ax_bot.set_xlim(-0.48, 2.48)
        ax_bot.set_ylim(-1.25 * q_max, 1.28 * q_max)
        ax_bot.set_xticks([0, 1, 2], q_labels, fontsize=7.8)
        ax_bot.set_ylabel("Four-tissue contrast (debiased UCell)")
        ax_bot.set_title(f"{title_base} · GSE164522 quartet Q", loc="left", fontsize=9.6, weight="bold")
        ax_bot.set_xlabel("Same four-tissue patients (PT, LM, NC, NL)")

    # Final column: B/Plasma::UPR falsification & B-vs-Plasma subtype decomposition reference
    upr_col = n_cols - 1
    upr_color = "#6b7280"
    ax_upr_top = axes[0, upr_col]
    ax_upr_top.axhline(0, color="#9ca3af", lw=0.85, ls="--", zorder=1)
    ax_upr_top.axvline(2.5, color="#cbd5e1", lw=1.0, ls=":", zorder=1)

    decomp_export = []
    d245 = pd.DataFrame()
    if m3_decomp_df is not None and not m3_decomp_df.empty:
        min_col = (
            "min_subtype_cells_per_tissue"
            if "min_subtype_cells_per_tissue" in m3_decomp_df.columns
            else "min_subtype_cells"
        )
        d245 = m3_decomp_df.loc[
            m3_decomp_df["program"].eq("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
            & m3_decomp_df[min_col].eq(1)
            & m3_decomp_df["status"].eq("complete_support")
        ].copy()
        for r in d245.itertuples(index=False):
            decomp_export.append(
                {
                    "dataset": "GSE245552",
                    "donor_uid": r.donor_uid,
                    "total_D_T": float(r.observed_total_D_T),
                    "composition": float(r.full_composition),
                    "within_subtype_state": float(r.full_within_subtype_state),
                }
            )

    d164 = pd.DataFrame()
    if m6_decomp_df is not None and not m6_decomp_df.empty:
        min_col164 = (
            "min_subtype_cells"
            if "min_subtype_cells" in m6_decomp_df.columns
            else "min_subtype_cells_per_tissue"
        )
        tot_col164 = (
            "observed_total_difference"
            if "observed_total_difference" in m6_decomp_df.columns
            else "observed_total_D_T"
        )
        st_col164 = (
            "full_within_type_state"
            if "full_within_type_state" in m6_decomp_df.columns
            else "full_within_subtype_state"
        )
        d164 = m6_decomp_df.loc[
            m6_decomp_df["program"].eq("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
            & m6_decomp_df[min_col164].eq(1)
            & m6_decomp_df["status"].eq("complete_support")
        ].copy()
        for r in d164.itertuples(index=False):
            decomp_export.append(
                {
                    "dataset": "GSE164522",
                    "donor_uid": r.donor_uid,
                    "total_D_T": float(getattr(r, tot_col164)),
                    "composition": float(r.full_composition),
                    "within_subtype_state": float(getattr(r, st_col164)),
                }
            )

    if m6_178_decomp_df is not None and not m6_178_decomp_df.empty:
        min_col178 = (
            "min_subtype_cells"
            if "min_subtype_cells" in m6_178_decomp_df.columns
            else "min_subtype_cells_per_tissue"
        )
        tot_col178 = (
            "observed_total_difference"
            if "observed_total_difference" in m6_178_decomp_df.columns
            else "observed_total_D_T"
        )
        st_col178 = (
            "full_within_type_state"
            if "full_within_type_state" in m6_178_decomp_df.columns
            else "full_within_subtype_state"
        )
        d178 = m6_178_decomp_df.loc[
            m6_178_decomp_df["program"].eq("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
            & m6_178_decomp_df[min_col178].eq(1)
            & m6_178_decomp_df["status"].eq("complete_support")
        ].copy()
        for r in d178.itertuples(index=False):
            decomp_export.append(
                {
                    "dataset": "GSE178318",
                    "donor_uid": r.donor_uid,
                    "total_D_T": float(getattr(r, tot_col178)),
                    "composition": float(r.full_composition),
                    "within_subtype_state": float(getattr(r, st_col178)),
                }
            )

    decomp_df_out = pd.DataFrame(decomp_export)
    if data_dir is not None and not decomp_df_out.empty:
        decomp_df_out.to_csv(data_dir / "upr_decomposition_reference.csv", index=False)

    upr_slots = [
        ("GSE245552", "total_D_T", "245\nTot", "#4b5563"),
        ("GSE245552", "composition", "245\nCmp", "#b45309"),
        ("GSE245552", "within_subtype_state", "245\nSt", "#2563eb"),
        ("GSE164522", "total_D_T", "164\nTot", "#4b5563"),
        ("GSE164522", "composition", "164\nCmp", "#b45309"),
        ("GSE164522", "within_subtype_state", "164\nSt", "#2563eb"),
    ]
    if not decomp_df_out.empty and "GSE178318" in set(decomp_df_out["dataset"]):
        ax_upr_top.axvline(5.5, color="#cbd5e1", lw=1.0, ls=":", zorder=1)
        ax_upr_top.axvspan(5.5, 8.5, color="#ecfdf5", alpha=0.65, zorder=0)
        upr_slots.extend(
            [
                ("GSE178318", "total_D_T", "178\nTot", "#4b5563"),
                ("GSE178318", "composition", "178\nCmp", "#b45309"),
                ("GSE178318", "within_subtype_state", "178\nSt", "#2563eb"),
            ]
        )
    upr_ticks = []
    for x, (ds, col_name, lbl, c_hex) in enumerate(upr_slots):
        vals = (
            decomp_df_out.loc[decomp_df_out["dataset"].eq(ds), col_name].dropna().to_numpy(dtype=float)
            if not decomp_df_out.empty
            else np.array([], dtype=float)
        )
        _scatter_with_mean(ax_upr_top, x, vals, c_hex, annotate_signs=False)
        upr_ticks.append(f"{lbl}\n({len(vals)})")
    upr_all = (
        decomp_df_out[["total_D_T", "composition", "within_subtype_state"]].to_numpy(dtype=float).ravel()
        if not decomp_df_out.empty
        else np.array([0.04])
    )
    upr_max = max(float(np.nanmax(np.abs(upr_all))), 0.03)
    ax_upr_top.set_xlim(-0.5, len(upr_slots) - 0.5)
    ax_upr_top.set_ylim(-1.22 * upr_max, 1.22 * upr_max)
    ax_upr_top.set_xticks(range(len(upr_slots)), upr_ticks, fontsize=6.6 if len(upr_slots) > 6 else 7.2)
    ax_upr_top.set_ylabel("Component of D_T (B vs Plasma)")
    ax_upr_top.set_title("B/Plasma · UPR (B vs Plasma decomposition; 203608 N/A)", loc="left", fontsize=9.1, weight="bold")
    ax_upr_top.set_xlabel("Symmetric decomposition (min ≥1/side; 178 N=4: 65.6% Cmp, 34.4% St)", fontsize=7.5)

    # Bottom right: B/Plasma::UPR GSE164522 Quartet Q
    ax_upr_bot = axes[1, upr_col]
    ax_upr_bot.axhline(0, color="#9ca3af", lw=0.85, ls="--", zorder=1)
    pp_upr = sub_patients.loc[
        sub_patients["lineage"].eq("B/Plasma")
        & sub_patients["program"].eq("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
    ]
    q_upr = pp_upr.loc[pp_upr["dataset"].eq("GSE164522") & pp_upr["I"].notna()].sort_values("donor_uid")
    if not q_upr.empty:
        for r in q_upr.itertuples(index=False):
            ax_upr_bot.plot(
                [0, 1, 2],
                [float(r.D_T_Q), float(r.D_N_Q), float(r.I)],
                color=upr_color,
                alpha=0.20,
                lw=0.85,
                zorder=2,
            )
    upr_q_ticks = []
    for x, (c_col, c_lbl) in enumerate(
        (("D_T_Q", "LM − PT\n(D_T_Q)"), ("D_N_Q", "NL − NC\n(D_N_Q)"), ("I", "Interaction\n(I = D_T − D_N)"))
    ):
        vals = q_upr[c_col].dropna().to_numpy(dtype=float) if not q_upr.empty else np.array([], dtype=float)
        _scatter_with_mean(ax_upr_bot, x, vals, upr_color)
        upr_q_ticks.append(f"{c_lbl}\nN = {len(vals)}")
    uq_all = (
        q_upr[["D_T_Q", "D_N_Q", "I"]].to_numpy(dtype=float).ravel()
        if not q_upr.empty
        else np.array([0.03])
    )
    uq_max = max(float(np.nanmax(np.abs(uq_all))), 0.03)
    ax_upr_bot.set_xlim(-0.48, 2.48)
    ax_upr_bot.set_ylim(-1.25 * uq_max, 1.28 * uq_max)
    ax_upr_bot.set_xticks([0, 1, 2], upr_q_ticks, fontsize=7.8)
    ax_upr_bot.set_ylabel("Four-tissue contrast (debiased UCell)")
    ax_upr_bot.set_title("B/Plasma · UPR · GSE164522 quartet Q", loc="left", fontsize=9.5, weight="bold")
    ax_upr_bot.set_xlabel("Same four-tissue patients (falsification check)")

    return save_or_verify(fig, "01_patient_programs", fig_dir)


def render_development_classifiers(
    metrics_df: pd.DataFrame,
    perm_null_df: pd.DataFrame | None,
    *,
    fig_dir: Path | None = FIG,
    data_dir: Path | None = DATA,
) -> dict[str, int]:
    """Build `02_development_classifiers` (`1 x 2`: single-lineage + All3Lineages with 200-perm null)."""
    base = metrics_df.loc[metrics_df["seed"].eq(20261121)].copy()
    if data_dir is not None:
        base.to_csv(data_dir / "development_classifier_metrics.csv", index=False)
        if perm_null_df is not None and not perm_null_df.empty:
            perm_null_df.to_csv(data_dir / "development_permutation_null.csv", index=False)

    fig, (ax_single, ax_joint) = plt.subplots(
        1, 2, figsize=(12.4, 6.1), gridspec_kw={"width_ratios": [1.28, 1.0]}, layout="constrained"
    )

    # Panel A: Single-lineage full-bag classifiers + 200-permutation empirical 95% null band
    present_single = [
        (m, lbl)
        for m, lbl in SINGLE_LINEAGE_METHODS
        if m in set(base.loc[base["lineage"].isin(["T/NK", "Mono/Mac/DC", "B/Plasma"]), "method"])
    ]
    if not present_single:
        raise ValueError("No single-lineage classifier rows found in patient_metrics.csv for seed 20261121")

    if perm_null_df is not None and not perm_null_df.empty:
        lo_env = float(perm_null_df["null_auroc_q025"].min())
        hi_env = float(perm_null_df["null_auroc_q975"].max())
        ax_single.axvspan(
            lo_env,
            hi_env,
            color="#e5e7eb",
            alpha=0.85,
            zorder=0,
            label=f"200-perm 95% null envelope [{lo_env:.2f}, {hi_env:.2f}]",
        )

    s_keys = [m for m, _ in present_single]
    s_labels = [lbl for _, lbl in present_single]
    for i, lineage in enumerate(["T/NK", "Mono/Mac/DC", "B/Plasma"]):
        sub = base.loc[base["lineage"].eq(lineage)].set_index("method").reindex(s_keys)
        y_pos = np.arange(len(s_keys)) + 0.20 * (i - 1)
        valid = sub["auroc"].notna().to_numpy()
        ax_single.scatter(
            sub.loc[valid, "auroc"],
            y_pos[valid],
            s=30,
            color=PALETTE[i],
            label=lineage,
            zorder=3,
        )

    ax_single.set_yticks(range(len(s_labels)), s_labels, fontsize=8.2)
    ax_single.invert_yaxis()
    ax_single.axvline(0.5, lw=0.9, ls="--", color="#6b7280", zorder=1)
    ax_single.set_xlim(0.0, 1.04)
    ax_single.set_xlabel("Held-out patient-site OOF AUROC (full clean bags)")
    ax_single.set_title(
        "A · GSE245552 single-lineage classifiers (fixed donor splits)",
        loc="left",
        weight="bold",
        fontsize=10.5,
    )
    ax_single.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.9, fontsize=7.8)

    # Panel B: All3Lineages cross-lineage joint specimen classifiers
    joint_df = base.loc[base["lineage"].eq("All3Lineages")].copy()
    present_joint = [
        (m, lbl) for m, lbl in JOINT_LINEAGE_METHODS if m in set(joint_df["method"])
    ]
    ax_joint.axvline(0.5, lw=0.9, ls="--", color="#6b7280", zorder=1)
    if present_joint:
        j_keys = [m for m, _ in present_joint]
        j_labels = [lbl for _, lbl in present_joint]
        j_sub = joint_df.set_index("method").reindex(j_keys)
        y_j = np.arange(len(j_keys))
        valid_j = j_sub["auroc"].notna().to_numpy()
        ax_joint.scatter(
            j_sub.loc[valid_j, "auroc"],
            y_j[valid_j],
            s=38,
            color=PALETTE[3],
            zorder=3,
            label="All3Lineages joint bag",
        )
        for idx_y, auc_v in zip(y_j[valid_j], j_sub.loc[valid_j, "auroc"]):
            ax_joint.text(
                float(auc_v) + 0.018,
                float(idx_y),
                f"{float(auc_v):.3f}",
                va="center",
                ha="left",
                fontsize=7.5,
                color="#1f2937",
            )
        n_donors_j = int(j_sub["patients"].max()) if "patients" in j_sub.columns else 16
        n_sites_j = int(j_sub["n_patient_sites"].max()) if "n_patient_sites" in j_sub.columns else 29
        ax_joint.set_yticks(range(len(j_labels)), j_labels, fontsize=8.2)
        ax_joint.invert_yaxis()
        ax_joint.set_title(
            f"B · All3Lineages joint specimen MIL ({n_donors_j} donors, {n_sites_j} sites)",
            loc="left",
            weight="bold",
            fontsize=10.5,
        )
    else:
        ax_joint.text(0.5, 0.5, "All3Lineages rows unavailable", ha="center", va="center")
        ax_joint.set_title("B · All3Lineages joint specimen MIL", loc="left", weight="bold", fontsize=10.5)

    ax_joint.set_xlim(0.0, 1.13)
    ax_joint.set_xlabel("Held-out patient-site OOF AUROC")
    return save_or_verify(fig, "02_development_classifiers", fig_dir)


def render_spatial_distances(
    distances_df: pd.DataFrame,
    support_df: pd.DataFrame,
    *,
    fig_dir: Path | None = FIG,
    data_dir: Path | None = DATA,
) -> dict[str, int]:
    """Build `03_spatial_distances` (`1 x 2`: 400->50 µm multi-scale decay + EPCAM+ nest stratification)."""
    if data_dir is not None:
        distances_df.to_csv(data_dir / "spatial_local_distances.csv", index=False)
        support_df.to_csv(data_dir / "spatial_sample_support.csv", index=False)

    cd8_map = {
        str(r.sample): int(r.n_cd3_cd8_primary) for r in support_df.itertuples(index=False)
    }
    my_pos_map = {
        str(r.sample): 100.0 * float(r.n_myeloid_primary_spp1) / max(int(r.n_myeloid_primary), 1)
        for r in support_df.itertuples(index=False)
    }

    fig, (ax_scale, ax_epcam) = plt.subplots(
        1, 2, figsize=(12.6, 4.9), gridspec_kw={"width_ratios": [1.05, 1.25]}, layout="constrained"
    )

    # Panel A: Multi-scale 400 -> 200 -> 100 -> 50 µm curves (dominance_disambiguated, all_myeloid, SPP1 >= 1)
    scales = [400, 200, 100, 50]
    x_pos = np.arange(len(scales), dtype=float)
    sub_scale = distances_df.loc[
        distances_df["branch"].eq("dominance_disambiguated")
        & distances_df["epcam_stratum"].eq("all_myeloid")
        & distances_df["spp1_min_transcripts"].eq(1)
    ].copy()

    section_order = sorted(
        sub_scale["sample"].unique(),
        key=lambda s: (0 if "AdjNormal" not in s else 1, s),
    )
    sec_colors = ["#175b75", "#b14925", "#577644", "#7c3aed"]
    ax_scale.axhline(0, color="#4b5563", lw=0.85, ls="--", zorder=1)

    for idx_s, sample_name in enumerate(section_order):
        rows_s = sub_scale.loc[sub_scale["sample"].eq(sample_name)].set_index("tile_um").reindex(scales)
        donor_short = str(rows_s["donor_uid"].dropna().iloc[0]).split(":")[-1]
        tissue = str(rows_s["tissue"].dropna().iloc[0])
        n_cd8 = cd8_map.get(sample_name, 0)
        t_min = int(rows_s["overlap_tiles"].min())
        t_max = int(rows_s["overlap_tiles"].max())
        lbl = f"{donor_short} · {tissue} (CD8 T={n_cd8:,}; {t_min}–{t_max} tiles)"
        col = sec_colors[idx_s % len(sec_colors)]
        jitter = (idx_s - 1.5) * 0.045

        y_vals = rows_s["delta_distance_um"].to_numpy(dtype=float)
        lo_vals = rows_s["null_lo"].to_numpy(dtype=float)
        hi_vals = rows_s["null_hi"].to_numpy(dtype=float)
        ls_style = "-" if tissue == "LM" else ":"
        marker = "o" if tissue == "LM" else "s"

        for k in range(len(scales)):
            if np.isfinite(lo_vals[k]) and np.isfinite(hi_vals[k]):
                ax_scale.plot(
                    [x_pos[k] + jitter, x_pos[k] + jitter],
                    [lo_vals[k], hi_vals[k]],
                    color=col,
                    alpha=0.30,
                    lw=3.2,
                    solid_capstyle="round",
                    zorder=2,
                )
        ax_scale.plot(
            x_pos + jitter,
            y_vals,
            color=col,
            ls=ls_style,
            lw=1.7,
            marker=marker,
            ms=5.5,
            label=lbl,
            zorder=4,
        )

    ax_scale.set_xticks(x_pos, [f"{t} µm" for t in scales])
    ax_scale.set_ylim(-44.0, 34.0)
    ax_scale.set_xlabel("Within-tile conditioning scale (coarse → fine)")
    ax_scale.set_ylabel("Nearest CD3/CD8 distance Δ (SPP1+ minus SPP1−, µm)")
    ax_scale.set_title(
        "A · Multi-scale decay & 95% permutation null bars",
        loc="left",
        weight="bold",
        fontsize=10.5,
    )
    ax_scale.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.92, fontsize=7.2)

    # Panel B: Tumor-nest (<15 µm to EPCAM+) vs stromal band (>=15 µm) at 200 µm and 100 µm across LM sections
    sub_ep = distances_df.loc[
        distances_df["branch"].eq("dominance_disambiguated")
        & distances_df["tissue"].eq("LM")
        & distances_df["spp1_min_transcripts"].eq(1)
        & distances_df["tile_um"].isin([200, 100])
        & distances_df["epcam_stratum"].isin(["near_epcam_lt15um", "stromal_epcam_ge15um"])
    ].copy()

    lm_samples = [s for s in section_order if "AdjNormal" not in s]
    ep_rows = []
    for s in lm_samples:
        for tile_w in (200, 100):
            for st_key, st_short in (
                ("near_epcam_lt15um", "<15 µm EPCAM+"),
                ("stromal_epcam_ge15um", "≥15 µm stroma"),
            ):
                hit = sub_ep.loc[
                    sub_ep["sample"].eq(s)
                    & sub_ep["tile_um"].eq(tile_w)
                    & sub_ep["epcam_stratum"].eq(st_key)
                ]
                if not hit.empty:
                    r0 = hit.iloc[0]
                    donor_short = str(r0["donor_uid"]).split(":")[-1]
                    ep_rows.append(
                        {
                            "sample": s,
                            "label": f"{donor_short} ({tile_w}µm) · {st_short}",
                            "stratum": st_key,
                            "delta": float(r0["delta_distance_um"]),
                            "null_lo": float(r0["null_lo"]),
                            "null_hi": float(r0["null_hi"]),
                            "n_tiles": int(r0["overlap_tiles"]),
                            "perm_p": float(r0["perm_p_two_sided"]),
                            "n_cd8": cd8_map.get(s, 0),
                            "spp1_pct": my_pos_map.get(s, np.nan),
                        }
                    )

    ep_df = pd.DataFrame(ep_rows)
    ax_epcam.axvline(0, color="#4b5563", lw=0.85, ls="--", zorder=1)
    for y_idx, r in enumerate(ep_df.itertuples(index=False)):
        ax_epcam.plot(
            [r.null_lo, r.null_hi],
            [y_idx, y_idx],
            color="#cbd5e1",
            lw=4.5,
            solid_capstyle="round",
            zorder=2,
            label="Within-tile 95% perm null" if y_idx == 0 else None,
        )
        is_nest = r.stratum == "near_epcam_lt15um"
        c_pt = "#b14925" if is_nest else "#175b75"
        lbl_pt = (
            "Near EPCAM+ tumor nest (<15 µm)"
            if (is_nest and y_idx == 0)
            else ("Stromal band (≥15 µm)" if (not is_nest and y_idx == 1) else None)
        )
        ax_epcam.scatter(
            r.delta,
            y_idx,
            s=40,
            color=c_pt,
            marker="D" if is_nest else "o",
            zorder=4,
            label=lbl_pt,
        )
        x_txt = max(r.delta, r.null_hi) + 0.6
        ax_epcam.text(
            x_txt,
            y_idx,
            f"{r.delta:+.2f} µm ({r.n_tiles} tiles, p={r.perm_p:.3f})",
            va="center",
            ha="left",
            fontsize=7.1,
            color="#1f2937",
        )
        if y_idx % 4 == 3 and y_idx < len(ep_df) - 1:
            ax_epcam.axhline(y_idx + 0.5, color="#e5e7eb", lw=0.8)

    if not ep_df.empty:
        max_x = float(np.nanmax(np.maximum(ep_df["delta"], ep_df["null_hi"])))
        min_x = float(np.nanmin(np.minimum(ep_df["delta"], ep_df["null_lo"])))
        ax_epcam.set_xlim(min_x - 2.0, max_x + 17.5)
        ax_epcam.set_yticks(range(len(ep_df)), ep_df["label"].tolist(), fontsize=7.7)
        ax_epcam.invert_yaxis()
    ax_epcam.set_xlabel("Nearest CD3/CD8 distance Δ (SPP1+ minus SPP1−, µm)")
    ax_epcam.set_title(
        "B · Tumor-nest vs stromal stratification (LM sections, SPP1 ≥ 1)",
        loc="left",
        weight="bold",
        fontsize=10.5,
    )
    ax_epcam.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=3,
        frameon=False,
        fontsize=7.5,
    )

    return save_or_verify(fig, "03_spatial_distances", fig_dir)


def _plot_external_cohort_row(
    row_axes: np.ndarray,
    v11_tab: pd.DataFrame,
    tf_tab: pd.DataFrame,
    v1_tab: pd.DataFrame,
    cohort_prefix: str,
    xlabel_text: str,
    *,
    include_legend_labels: bool,
) -> None:
    m_keys = [m for m, _ in EXTERNAL_METHODS]
    m_labels = [lbl for _, lbl in EXTERNAL_METHODS]
    for idx_l, (ax, lineage, color) in enumerate(
        zip(row_axes, ["T/NK", "Mono/Mac/DC", "B/Plasma"], PALETTE[:3])
    ):
        r_v11 = v11_tab.loc[v11_tab["lineage"].eq(lineage)].set_index("method").reindex(m_keys)
        r_tf = (
            tf_tab.loc[tf_tab["lineage"].eq(lineage)].set_index("method").reindex(m_keys)
            if not tf_tab.empty
            else None
        )
        r_v1 = (
            v1_tab.loc[v1_tab["lineage"].eq(lineage)].set_index("method").reindex(m_keys)
            if not v1_tab.empty
            else None
        )

        for i, m in enumerate(m_keys):
            row11 = r_v11.loc[m]
            y11 = i - 0.15 if (r_tf is not None or r_v1 is not None) else float(i)
            if pd.notna(row11.get("auroc", np.nan)):
                ax.plot(
                    [float(row11["auroc_lower95"]), float(row11["auroc_upper95"])],
                    [y11, y11],
                    color=color,
                    lw=2.0,
                    alpha=0.85,
                    zorder=3,
                )
                ax.scatter(
                    float(row11["auroc"]),
                    y11,
                    color=color,
                    s=36,
                    zorder=4,
                    label="V1.1 (full clean bags, 15/15)" if (include_legend_labels and idx_l == 0 and i == 0) else None,
                )
            else:
                ax.text(0.5, y11, "not estimable", ha="center", va="center", fontsize=7.5)

            if r_tf is not None and m in r_tf.index and pd.notna(r_tf.loc[m].get("auroc", np.nan)):
                row_tf = r_tf.loc[m]
                ytf = i + 0.03
                ax.plot(
                    [float(row_tf["auroc_lower95"]), float(row_tf["auroc_upper95"])],
                    [ytf, ytf],
                    color="#6b7280",
                    lw=1.0,
                    ls=":",
                    alpha=0.75,
                    zorder=2,
                )
                ax.scatter(
                    float(row_tf["auroc"]),
                    ytf,
                    facecolors="white",
                    edgecolors="#4b5563",
                    marker="D",
                    s=22,
                    lw=1.0,
                    zorder=3,
                    label="V1_F1_techfix (15/15)" if (include_legend_labels and idx_l == 0 and i == 0) else None,
                )

            if r_v1 is not None and m in r_v1.index:
                row_v1 = r_v1.loc[m]
                yv1 = i + 0.20
                if pd.notna(row_v1.get("auroc", np.nan)):
                    ax.scatter(
                        float(row_v1["auroc"]),
                        yv1,
                        color="#9ca3af",
                        marker="x",
                        s=24,
                        lw=1.1,
                        zorder=2,
                        label="V1_F1 baseline (11/15, 4 abstained)" if (include_legend_labels and idx_l == 0 and i == 0) else None,
                    )
                else:
                    ax.text(
                        0.50,
                        yv1,
                        "V1 abstained (<95% raw HVG)",
                        ha="center",
                        va="center",
                        fontsize=6.5,
                        color="#b91c1c",
                    )

        ax.axvline(0.5, ls="--", color="#9ca3af", lw=0.85, zorder=1)
        ax.set_xlim(-0.03, 1.04)
        n_don = int(r_v11["n_donors"].max()) if r_v11["n_donors"].notna().any() else 0
        n_sit = int(r_v11["n_sites"].max()) if r_v11["n_sites"].notna().any() else 0
        ax.set_title(
            f"{cohort_prefix} · {lineage} ({n_don}d, {n_sit}s)",
            loc="left",
            weight="bold",
            fontsize=9.0,
        )
        ax.set_xlabel(xlabel_text)

    row_axes[0].set_yticks(range(len(m_labels)), m_labels, fontsize=8.3)


def _collect_three_version_rows(
    cohort: str,
    v1_tab: pd.DataFrame,
    tf_tab: pd.DataFrame,
    v11_tab: pd.DataFrame,
) -> list[dict[str, object]]:
    m_keys = [m for m, _ in EXTERNAL_METHODS]
    out: list[dict[str, object]] = []
    for ver_name, tab in (
        ("CRLM_V1_245_F1", v1_tab),
        ("CRLM_V1_245_F1_techfix", tf_tab),
        ("CRLM_V1.1_245_F1", v11_tab),
    ):
        for lineage in ("T/NK", "Mono/Mac/DC", "B/Plasma"):
            sub = (
                tab.loc[tab["lineage"].eq(lineage)].set_index("method").reindex(m_keys)
                if not tab.empty
                else pd.DataFrame(index=m_keys)
            )
            for m in m_keys:
                row = sub.loc[m] if m in sub.index else pd.Series(dtype=object)
                auc = row.get("auroc", np.nan)
                status = str(row.get("status", "")).strip()
                if status in ("", "nan", "None"):
                    if pd.notna(auc):
                        status = "predicted"
                    elif m.startswith("PCA_") and pd.notna(row.get("n_sites", np.nan)) and int(row["n_sites"]) == 0:
                        status = "abstained_lt95pct_raw_hvg"
                    else:
                        status = "not_estimable"
                out.append(
                    {
                        "cohort": cohort,
                        "version": ver_name,
                        "lineage": lineage,
                        "method": m,
                        "status": status,
                        "n_donors": int(row["n_donors"]) if pd.notna(row.get("n_donors", np.nan)) else 0,
                        "n_sites": int(row["n_sites"]) if pd.notna(row.get("n_sites", np.nan)) else 0,
                        "paired_n": int(row["paired_n"]) if pd.notna(row.get("paired_n", np.nan)) else 0,
                        "auroc": float(auc) if pd.notna(auc) else np.nan,
                        "auroc_lower95": float(row["auroc_lower95"]) if pd.notna(row.get("auroc_lower95", np.nan)) else np.nan,
                        "auroc_upper95": float(row["auroc_upper95"]) if pd.notna(row.get("auroc_upper95", np.nan)) else np.nan,
                        "average_precision": float(row["average_precision"]) if pd.notna(row.get("average_precision", np.nan)) else np.nan,
                        "balanced_accuracy": float(row["balanced_accuracy"]) if pd.notna(row.get("balanced_accuracy", np.nan)) else np.nan,
                        "fraction_LM_above_PT": float(row["fraction_LM_above_PT"]) if pd.notna(row.get("fraction_LM_above_PT", np.nan)) else np.nan,
                    }
                )
    return out


def render_external_classifiers(
    v11_metrics_df: pd.DataFrame,
    techfix_metrics_df: pd.DataFrame | None = None,
    v1_archive_metrics_df: pd.DataFrame | None = None,
    v11_178_metrics_df: pd.DataFrame | None = None,
    three_ver_178_df: pd.DataFrame | None = None,
    *,
    fig_dir: Path | None = FIG,
    data_dir: Path | None = DATA,
) -> dict[str, int]:
    """Build `04_external_classifiers` (`2 x 3`: GSE164522 + prospective GSE178318 15/15 transfer)."""
    if data_dir is not None:
        v11_metrics_df.to_csv(data_dir / "external_classifier_metrics.csv", index=False)
        if v11_178_metrics_df is not None and not v11_178_metrics_df.empty:
            v11_178_metrics_df.to_csv(data_dir / "external_gse178318_classifier_metrics.csv", index=False)

    v11_tab = v11_metrics_df.loc[v11_metrics_df["scope"].eq("common_predicting_methods")].copy()
    tf_tab = pd.DataFrame()
    v1_tab = pd.DataFrame()
    if techfix_metrics_df is not None and not techfix_metrics_df.empty:
        tf_tab = techfix_metrics_df.loc[
            techfix_metrics_df["step"].eq("step_b_full_techfix")
            & techfix_metrics_df["scope"].eq("method_available")
        ].copy()
        if v1_archive_metrics_df is None or v1_archive_metrics_df.empty:
            v1_tab = techfix_metrics_df.loc[
                techfix_metrics_df["step"].eq("m6_original_published")
                & techfix_metrics_df["scope"].eq("method_available")
            ].copy()
    if v1_tab.empty and v1_archive_metrics_df is not None and not v1_archive_metrics_df.empty:
        v1_tab = v1_archive_metrics_df.loc[
            v1_archive_metrics_df["scope"].eq("method_available")
        ].copy()

    plotted_3ver_rows = _collect_three_version_rows("GSE164522", v1_tab, tf_tab, v11_tab)

    has_178 = v11_178_metrics_df is not None and not v11_178_metrics_df.empty
    n_rows = 2 if has_178 else 1
    fig, axes = plt.subplots(
        n_rows,
        3,
        figsize=(12.4, 7.2 if has_178 else 4.2),
        sharey=True,
        layout="constrained",
    )
    axes_2d = np.atleast_2d(axes)

    _plot_external_cohort_row(
        axes_2d[0],
        v11_tab,
        tf_tab,
        v1_tab,
        "A · 164 (Pri)",
        "GSE164522 patient-site AUROC",
        include_legend_labels=True,
    )

    if has_178:
        v11_178_tab = v11_178_metrics_df.loc[v11_178_metrics_df["scope"].eq("method_available")].copy()
        tf_178_tab = pd.DataFrame()
        v1_178_tab = pd.DataFrame()
        if three_ver_178_df is not None and not three_ver_178_df.empty:
            tf_178_tab = three_ver_178_df.loc[three_ver_178_df["version"].eq("CRLM_V1_245_F1_techfix")].copy()
            v1_178_tab = three_ver_178_df.loc[three_ver_178_df["version"].eq("CRLM_V1_245_F1")].copy()
        plotted_3ver_rows.extend(_collect_three_version_rows("GSE178318", v1_178_tab, tf_178_tab, v11_178_tab))
        _plot_external_cohort_row(
            axes_2d[1],
            v11_178_tab,
            tf_178_tab,
            v1_178_tab,
            "B · 178 (Zero-Exp)",
            "GSE178318 prospective patient-site AUROC",
            include_legend_labels=False,
        )
    if data_dir is not None:
        pd.DataFrame(plotted_3ver_rows).to_csv(
            data_dir / "three_version_classifier_comparison.csv", index=False
        )
    axes_2d[0, 0].invert_yaxis()

    handles, labels = axes_2d[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.05 if has_178 else -0.08),
            ncol=len(handles),
            frameon=False,
            fontsize=8.2,
        )
    fig.suptitle(
        "Dual-cohort external transfer (GSE164522 + prospective GSE178318) · V1.1 100% HVG coverage (15/15 active)",
        weight="bold",
        fontsize=11.0,
    )
    return save_or_verify(fig, "04_external_classifiers", fig_dir)


def build_all_from_disk(root_dir: Path = ROOT) -> dict[str, dict[str, int]]:
    """Read stage tables from `artifacts/{m3,m4,m5,m6,m7}` and write `artifacts/m8/{figures,source_data}`."""
    apply_rcparams()
    art = root_dir / "artifacts"
    out_m8 = art / "m8"
    fig_dir = out_m8 / "figures"
    data_dir = out_m8 / "source_data"
    fig_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    patients = pd.concat(
        [
            pd.read_csv(art / "m3/patient_effects/patient_contrasts.csv"),
            pd.read_csv(art / "m6/patient_contrasts.csv"),
        ],
        ignore_index=True,
    )
    cand_path = art / "m5/frozen_candidate_registry.csv"
    candidates = pd.read_csv(cand_path) if cand_path.is_file() else None
    m3_decomp_path = art / "m3/composition_results/patient_subtype_decomposition.csv"
    m3_decomp = pd.read_csv(m3_decomp_path) if m3_decomp_path.is_file() else None
    m6_decomp_path = art / "m6/GSE164522/patient_effects/patient_B_plasma_decomposition.csv"
    m6_decomp = pd.read_csv(m6_decomp_path) if m6_decomp_path.is_file() else None
    m6_178_decomp_path = art / "m6/GSE178318/patient_effects/patient_B_plasma_decomposition.csv"
    m6_178_decomp = pd.read_csv(m6_178_decomp_path) if m6_178_decomp_path.is_file() else None

    m4_metrics = pd.read_csv(art / "m4/patient_metrics.csv")
    perm_path = art / "m4/PCA/permutation_null_summary.csv"
    perm_null = pd.read_csv(perm_path) if perm_path.is_file() else None

    spatial_dist = pd.read_csv(art / "m7/spatial_measurement_audit/local_distances.csv")
    spatial_sup = pd.read_csv(art / "m7/spatial_measurement_audit/sample_support.csv")

    ext_v11 = pd.read_csv(art / "m6/GSE164522/model_predictions/patient_metrics.csv")
    ext_tf = pd.read_csv(TECHFIX_EXT_METRICS) if TECHFIX_EXT_METRICS.is_file() else None
    ext_v1 = pd.read_csv(ARCHIVE_V1_EXT_METRICS) if ARCHIVE_V1_EXT_METRICS.is_file() else None
    ext_178_path = art / "m6/GSE178318/model_predictions/patient_metrics.csv"
    ext_178 = pd.read_csv(ext_178_path) if ext_178_path.is_file() else None
    ext_178_3v_path = art / "m6/GSE178318/model_predictions/three_version_gse178318_classifiers.csv"
    ext_178_3v = pd.read_csv(ext_178_3v_path) if ext_178_3v_path.is_file() else None

    return {
        "01_patient_programs": render_patient_programs(
            patients,
            candidates,
            m3_decomp,
            m6_decomp,
            m6_178_decomp,
            fig_dir=fig_dir,
            data_dir=data_dir,
        ),
        "02_development_classifiers": render_development_classifiers(
            m4_metrics, perm_null, fig_dir=fig_dir, data_dir=data_dir
        ),
        "03_spatial_distances": render_spatial_distances(
            spatial_dist, spatial_sup, fig_dir=fig_dir, data_dir=data_dir
        ),
        "04_external_classifiers": render_external_classifiers(
            ext_v11,
            ext_tf,
            ext_v1,
            ext_178,
            ext_178_3v,
            fig_dir=fig_dir,
            data_dir=data_dir,
        ),
    }


def self_check() -> None:
    """Render all 4 V1.1 figures purely in memory (`io.BytesIO`, zero disk writes)."""
    apply_rcparams()
    rng = np.random.default_rng(20261121)

    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "C1",
                "candidate_tier": "primary",
                "lineage": "Mono/Mac/DC",
                "program": "Core_Lipid_SAM_R",
                "expected_direction_DT": "positive",
            },
            {
                "candidate_id": "C2",
                "candidate_tier": "primary",
                "lineage": "T/NK",
                "program": "B07_Th17_TWEAK",
                "expected_direction_DT": "negative",
            },
        ]
    )

    pat_rows = []
    for ds, n_d, has_q in (("GSE245552", 11, False), ("GSE164522", 8, True), ("GSE203608", 4, False)):
        for i in range(n_d):
            for lin, prog, mu in (
                ("Mono/Mac/DC", "Core_Lipid_SAM_R", 0.12),
                ("T/NK", "B07_Th17_TWEAK", -0.03),
                ("B/Plasma", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", -0.01),
            ):
                dt = float(rng.normal(mu, 0.02))
                dn = float(rng.normal(0.01, 0.015)) if (has_q and i < 5) else np.nan
                pat_rows.append(
                    {
                        "dataset": ds,
                        "donor_uid": f"{ds}:D{i:02d}",
                        "lineage": lin,
                        "program": prog,
                        "D_T": dt,
                        "D_N": dn,
                        "D_T_Q": dt if np.isfinite(dn) else np.nan,
                        "D_N_Q": dn,
                        "I": (dt - dn) if np.isfinite(dn) else np.nan,
                    }
                )
    patients_syn = pd.DataFrame(pat_rows)

    m3_dec = pd.DataFrame(
        [
            {
                "donor_uid": f"GSE245552:D{i:02d}",
                "lineage": "B/Plasma",
                "program": "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                "min_subtype_cells_per_tissue": 1,
                "status": "complete_support",
                "observed_total_D_T": 0.008 + 0.002 * i,
                "full_composition": 0.009 + 0.001 * i,
                "full_within_subtype_state": -0.001 + 0.001 * i,
            }
            for i in range(8)
        ]
    )
    m6_dec = pd.DataFrame(
        [
            {
                "donor_uid": f"GSE164522:D{i:02d}",
                "lineage": "B/Plasma",
                "program": "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                "min_subtype_cells": 1,
                "status": "complete_support",
                "observed_total_difference": -0.015 + 0.002 * i,
                "full_composition": -0.013 + 0.001 * i,
                "full_within_type_state": -0.002 + 0.001 * i,
            }
            for i in range(8)
        ]
    )

    m4_rows = []
    for lin in ("T/NK", "Mono/Mac/DC", "B/Plasma"):
        for m, _ in SINGLE_LINEAGE_METHODS:
            m4_rows.append(
                {
                    "method": m,
                    "lineage": lin,
                    "seed": 20261121,
                    "patients": 15,
                    "n_patient_sites": 27,
                    "auroc": float(np.clip(rng.normal(0.84, 0.08), 0.45, 1.0)),
                }
            )
    for m, _ in JOINT_LINEAGE_METHODS:
        m4_rows.append(
            {
                "method": m,
                "lineage": "All3Lineages",
                "seed": 20261121,
                "patients": 16,
                "n_patient_sites": 29,
                "auroc": float(np.clip(rng.normal(0.91, 0.05), 0.5, 1.0)),
            }
        )
    perm_syn = pd.DataFrame(
        [
            {"lineage": lin, "null_auroc_q025": 0.24, "null_auroc_q975": 0.75}
            for lin in ("T/NK", "Mono/Mac/DC", "B/Plasma")
        ]
    )

    # Read real already-verified M7 & Techfix CSVs (read-only) to test exact schema fidelity
    spatial_dist = pd.read_csv(A / "m7/spatial_measurement_audit/local_distances.csv")
    spatial_sup = pd.read_csv(A / "m7/spatial_measurement_audit/sample_support.csv")
    ext_tf = pd.read_csv(TECHFIX_EXT_METRICS)
    ext_v11_syn = (
        ext_tf.loc[
            ext_tf["step"].eq("step_b_full_techfix")
            & ext_tf["scope"].eq("common_predicting_methods")
        ]
        .drop(columns=["step"])
        .copy()
    )

    s1 = render_patient_programs(
        patients_syn, candidates, m3_dec, m6_dec, fig_dir=None, data_dir=None
    )
    s2 = render_development_classifiers(
        pd.DataFrame(m4_rows), perm_syn, fig_dir=None, data_dir=None
    )
    s3 = render_spatial_distances(spatial_dist, spatial_sup, fig_dir=None, data_dir=None)
    s4 = render_external_classifiers(ext_v11_syn, ext_tf, None, fig_dir=None, data_dir=None)

    print(
        f"build_paper_figures --self-check PASSED (zero disk writes): "
        f"01={s1['png']}B, 02={s2['png']}B, 03={s3['png']}B, 04={s4['png']}B"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="Render all 4 V1.1 figures purely into io.BytesIO buffers (zero disk writes).",
    )
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        sizes = build_all_from_disk(ROOT)
        for name, sz in sizes.items():
            print(f"{name}: png={sz['png']}B, pdf={sz['pdf']}B")


if __name__ == "__main__":
    main()
