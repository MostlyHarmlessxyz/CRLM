"""Assemble V1.1 external validation results, apply family-wise Holm adjustment, and build the
three-version side-by-side evidence table (`V1_F1` vs `V1_F1_techfix` vs `V1.1`)."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "starter"))

from adjust_registered_tests import adjust_registered, run as run_adjust

ARCHIVE_V1_M6 = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m6/registered_tests/registered_test_results.csv"
TECHFIX_M6 = ROOT / "artifacts/techfix_v1_f1/programs/registered_test_results.csv"
ARCHIVE_V1_CLS = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m6/GSE164522/model_predictions/patient_metrics.csv"
TECHFIX_CLS = ROOT / "artifacts/techfix_v1_f1/models/techfix_external_classifier_metrics.csv"
V11_CLS = ROOT / "artifacts/m6/GSE164522/model_predictions/patient_metrics.csv"


def _sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _fmt_cell(row: pd.Series | None) -> str:
    if row is None:
        return "未登记该版本主检验"
    raw_p = row.get("raw_p", row.get("sign_p_two_sided", np.nan))
    if pd.isna(raw_p):
        n_pat = row.get("n_patients", np.nan)
        if pd.notna(n_pat) and int(n_pat) == 0:
            return "N=0 (覆盖不足)"
        return "NA (条款待核, Holm p=1.0)"
    mean_val = float(row["mean_difference"])
    n_pos = int(row["n_positive"])
    n_neg = int(row["n_negative"])
    n_pat = int(row["n_patients"])
    holm_p = row.get("holm_adjusted_p", np.nan)
    if pd.notna(holm_p):
        return f"{mean_val:+.5f} ({n_pos}/{n_neg}, N={n_pat}; p={float(raw_p):.4f}, p_H={float(holm_p):.4f})"
    return f"{mean_val:+.5f} ({n_pos}/{n_neg}, N={n_pat}; p={float(raw_p):.4f}, 探索性对照)"


def build_three_version_comparison(
    v1_df: pd.DataFrame,
    techfix_df: pd.DataFrame,
    v11_reg_df: pd.DataFrame,
    v11_expl_df: pd.DataFrame,
) -> tuple[pd.DataFrame, str]:
    """Build structured CSV & Markdown side-by-side table across V1_F1, V1_F1_techfix, and V1.1."""
    key_cols = ["cohort", "lineage", "program", "contrast"]

    def index_by_key(df: pd.DataFrame) -> dict[tuple[str, str, str, str], pd.Series]:
        out: dict[tuple[str, str, str, str], pd.Series] = {}
        if df.empty:
            return out
        cohort_col = "cohort" if "cohort" in df.columns else "dataset"
        for _, r in df.iterrows():
            k = (str(r[cohort_col]), str(r["lineage"]), str(r["program"]), str(r["contrast"]))
            if k not in out:
                out[k] = r
        return out

    m_v1 = index_by_key(v1_df)
    m_tf = index_by_key(techfix_df)
    m_v11_reg = index_by_key(v11_reg_df)
    m_v11_exp = index_by_key(v11_expl_df)

    # Also match B19_S4_Mac_SPP1_clean exploratory row to historical B19_S4_Mac_SPP1 when exact name differs
    ordered_keys: list[tuple[str, str, str, str]] = []
    for src in (v11_reg_df, v1_df):
        if src.empty:
            continue
        cohort_col = "cohort" if "cohort" in src.columns else "dataset"
        for _, r in src.iterrows():
            k = (str(r[cohort_col]), str(r["lineage"]), str(r["program"]), str(r["contrast"]))
            if k not in ordered_keys:
                ordered_keys.append(k)

    rows = []
    md_lines = [
        "# 三版本外部队列预登记检验并列对照表（`V1_F1` vs `V1_F1_techfix` vs `V1.1`）",
        "",
        "| 队列 | 对比 | 谱系 :: 程序 | `CRLM_V1_245_F1`（原版基线） | `CRLM_V1_245_F1_techfix`（纯技术勘误） | `CRLM_V1.1_245_F1`（方法学重构去偏） |",
        "|---|---|---|---|---|---|",
    ]

    for k in ordered_keys:
        cohort, lineage, prog, contrast = k
        r0 = m_v1.get(k)
        rtf = m_tf.get(k)
        r11 = m_v11_reg.get(k)
        v11_tier = "registered"
        if r11 is None:
            r11 = m_v11_exp.get(k)
            if r11 is None and prog == "B19_S4_Mac_SPP1":
                r11 = m_v11_reg.get((cohort, lineage, "B19_S4_Mac_SPP1_clean", contrast)) or m_v11_exp.get(
                    (cohort, lineage, "B19_S4_Mac_SPP1_clean", contrast)
                )
            v11_tier = "exploratory_sensitivity" if r11 is not None else "not_registered"

        c0_txt = _fmt_cell(r0)
        ctf_txt = _fmt_cell(rtf)
        c11_txt = _fmt_cell(r11)
        md_lines.append(f"| {cohort} | `{contrast}` | `{lineage}::{prog}` | {c0_txt} | {ctf_txt} | {c11_txt} |")

        rows.append(
            {
                "cohort": cohort,
                "contrast": contrast,
                "lineage": lineage,
                "program": prog,
                "v1_f1_summary": c0_txt,
                "v1_f1_mean": float(r0["mean_difference"]) if (r0 is not None and pd.notna(r0.get("mean_difference"))) else np.nan,
                "v1_f1_holm_p": float(r0["holm_adjusted_p"]) if (r0 is not None and pd.notna(r0.get("holm_adjusted_p"))) else np.nan,
                "techfix_summary": ctf_txt,
                "techfix_mean": float(rtf["mean_difference"]) if (rtf is not None and pd.notna(rtf.get("mean_difference"))) else np.nan,
                "techfix_holm_p": float(rtf["holm_adjusted_p"]) if (rtf is not None and pd.notna(rtf.get("holm_adjusted_p"))) else np.nan,
                "v1_1_tier": v11_tier,
                "v1_1_summary": c11_txt,
                "v1_1_mean": float(r11["mean_difference"]) if (r11 is not None and pd.notna(r11.get("mean_difference"))) else np.nan,
                "v1_1_raw_p": float(r11.get("raw_p", r11.get("sign_p_two_sided", np.nan))) if r11 is not None else np.nan,
                "v1_1_holm_p": float(r11["holm_adjusted_p"]) if (r11 is not None and "holm_adjusted_p" in r11 and pd.notna(r11["holm_adjusted_p"])) else np.nan,
            }
        )

    return pd.DataFrame(rows), "\n".join(md_lines) + "\n"


def _fmt_cls(row: pd.Series | None) -> str:
    if row is None or pd.isna(row.get("auroc", np.nan)):
        return "拒判 (<95% 原始符号 HVG)"
    auc = float(row["auroc"])
    lo = float(row["auroc_lower95"])
    hi = float(row["auroc_upper95"])
    n_p = int(row.get("paired_n", 0))
    frac = row.get("fraction_LM_above_PT", np.nan)
    pos_pairs = f"{int(round(float(frac) * n_p))}/{n_p}" if (n_p > 0 and pd.notna(frac)) else f"N={int(row['n_donors'])}"
    return f"{auc:.3f} [{lo:.3f}, {hi:.3f}] ({pos_pairs})"


def build_classifier_three_version_comparison(
    v1_cls: pd.DataFrame,
    tf_cls: pd.DataFrame,
    v11_cls: pd.DataFrame,
) -> tuple[pd.DataFrame, str]:
    """Build structured CSV & Markdown section comparing GSE164522 15 external transfer classifiers."""
    if v11_cls.empty:
        return pd.DataFrame(), ""
    v1_sub = v1_cls.loc[v1_cls["scope"].eq("method_available")] if not v1_cls.empty else pd.DataFrame()
    tf_sub = (
        tf_cls.loc[tf_cls["step"].eq("step_b_full_techfix") & tf_cls["scope"].eq("method_available")]
        if not tf_cls.empty
        else pd.DataFrame()
    )
    v11_sub = v11_cls.loc[v11_cls["scope"].eq("method_available")]

    def lookup(df: pd.DataFrame) -> dict[tuple[str, str], pd.Series]:
        return {(str(r["lineage"]), str(r["method"])): r for _, r in df.iterrows()}

    m0, mtf, m11 = lookup(v1_sub), lookup(tf_sub), lookup(v11_sub)
    methods = ("program_elastic", "PCA_mean", "PCA_attention", "scGPT_mean", "scGPT_attention")
    lineages = ("T/NK", "Mono/Mac/DC", "B/Plasma")

    rows = []
    md = [
        "",
        "## 二、GSE164522 外部冻结迁移分类器三版本对照（15 项模型，患者组织级 AUROC [95% CI] 与 LM>PT 配对数）",
        "",
        "| 谱系 | 冻结分类头 | `CRLM_V1_245_F1`（11/15 激活，4 拒判） | `CRLM_V1_245_F1_techfix`（15/15 激活） | `CRLM_V1.1_245_F1`（全袋洁净，15/15 激活） |",
        "|---|---|---|---|---|",
    ]
    for lin in lineages:
        for meth in methods:
            r0 = m0.get((lin, meth))
            rtf = mtf.get((lin, meth))
            r11 = m11.get((lin, meth))
            s0, stf, s11 = _fmt_cls(r0), _fmt_cls(rtf), _fmt_cls(r11)
            md.append(f"| `{lin}` | `{meth}` | {s0} | {stf} | {s11} |")
            rows.append(
                {
                    "lineage": lin,
                    "method": meth,
                    "v1_f1_auroc": float(r0["auroc"]) if (r0 is not None and pd.notna(r0.get("auroc"))) else np.nan,
                    "v1_f1_summary": s0,
                    "techfix_auroc": float(rtf["auroc"]) if (rtf is not None and pd.notna(rtf.get("auroc"))) else np.nan,
                    "techfix_summary": stf,
                    "v1_1_auroc": float(r11["auroc"]) if (r11 is not None and pd.notna(r11.get("auroc"))) else np.nan,
                    "v1_1_lower95": float(r11["auroc_lower95"]) if (r11 is not None and pd.notna(r11.get("auroc_lower95"))) else np.nan,
                    "v1_1_upper95": float(r11["auroc_upper95"]) if (r11 is not None and pd.notna(r11.get("auroc_upper95"))) else np.nan,
                    "v1_1_summary": s11,
                }
            )
    return pd.DataFrame(rows), "\n".join(md) + "\n"


def main() -> None:
    out = ROOT / "artifacts/m6"
    plan_path = ROOT / "artifacts/m5/registered_tests.csv"
    plan_df = pd.read_csv(plan_path)
    n_families = plan_df["family"].nunique()
    if len(plan_df) % n_families != 0:
        raise ValueError("Expected equal number of primary candidates across registered families")
    n_per_family = len(plan_df) // n_families

    tables = []
    patient_tables = []
    sens_tables = []
    focal_tables = []

    for accession in ("GSE164522", "GSE203608", "GSE178318"):
        folder = out / accession / "patient_effects"
        if accession == "GSE178318" and not (folder / "effects_summary.csv").exists():
            continue
        summary = pd.read_csv(folder / "effects_summary.csv")
        patients = pd.read_csv(folder / "patient_contrasts.csv")
        keys = ["dataset", "donor_uid", "lineage", "program"]
        if patients.duplicated(keys).any() or set(patients.dataset) != {accession}:
            raise ValueError(f"Invalid external patient identity for {accession}")
        q = patients[patients.I.notna()]
        if not q.empty and not (
            q[["A_PT", "A_LM", "A_NC", "A_NL"]].notna().all().all()
            and np.allclose(q.I, q.D_T_Q - q.D_N_Q)
            and np.allclose(q.D_T_Q, q.D_T)
            and np.allclose(q.D_N_Q, q.D_N)
        ):
            raise ValueError(f"Four-tissue contrasts differ from their common patients in {accession}")
        tables.append(summary)
        patient_tables.append(patients)

        sens_path = folder / "sensitivity_summary.csv"
        if sens_path.is_file():
            sens_tables.append(pd.read_csv(sens_path))
        focal_path = folder / "focal_subtype_effects_summary.csv"
        if focal_path.is_file():
            focal_tables.append(pd.read_csv(focal_path))

    summary_path = out / "effects_summary.csv"
    pd.concat(tables, ignore_index=True).to_csv(summary_path, index=False)
    pd.concat(patient_tables, ignore_index=True).to_csv(out / "patient_contrasts.csv", index=False)
    if sens_tables:
        sens_all = pd.concat(sens_tables, ignore_index=True)
        sens_all.to_csv(out / "sensitivity_summary_all_cohorts.csv", index=False)
        rule_sub = sens_all.loc[sens_all["analysis"].astype(str).str.startswith("rule_annotation_m1")]
        if not rule_sub.empty:
            rule_sub.to_csv(out / "uniform_annotation_effects_summary.csv", index=False)
    if focal_tables:
        pd.concat(focal_tables, ignore_index=True).to_csv(
            out / "focal_subtype_effects_all_cohorts.csv", index=False
        )

    run_adjust(plan_path, summary_path, out / "registered_tests")
    adjusted = pd.read_csv(out / "registered_tests/registered_test_results.csv")
    exploratory_path = out / "registered_tests/unplanned_exploratory_results.csv"
    exploratory = pd.read_csv(exploratory_path) if exploratory_path.is_file() else pd.DataFrame()

    assert len(adjusted) == len(plan_df) and adjusted.test_id.is_unique
    assert adjusted.groupby("family").size().eq(n_per_family).all()

    for accession in ("GSE203608", "GSE178318"):
        folder = out / accession / "registered_tests"
        folder.mkdir(parents=True, exist_ok=True)
        adjusted.loc[adjusted.cohort.eq(accession)].to_csv(folder / "registered_results.csv", index=False)

    v1_df = pd.read_csv(ARCHIVE_V1_M6) if ARCHIVE_V1_M6.is_file() else pd.DataFrame()
    tf_df = pd.read_csv(TECHFIX_M6) if TECHFIX_M6.is_file() else pd.DataFrame()
    comp_df, comp_md = build_three_version_comparison(v1_df, tf_df, adjusted, exploratory)

    v1_cls = pd.read_csv(ARCHIVE_V1_CLS) if ARCHIVE_V1_CLS.is_file() else pd.DataFrame()
    tf_cls = pd.read_csv(TECHFIX_CLS) if TECHFIX_CLS.is_file() else pd.DataFrame()
    v11_cls = pd.read_csv(V11_CLS) if V11_CLS.is_file() else pd.DataFrame()
    cls_df, cls_md = build_classifier_three_version_comparison(v1_cls, tf_cls, v11_cls)
    if not cls_df.empty:
        cls_df.to_csv(out / "three_version_classifier_comparison.csv", index=False)
        comp_md = comp_md + cls_md

    comp_df.to_csv(out / "three_version_evidence_comparison.csv", index=False)
    (out / "three_version_evidence_comparison.md").write_text(comp_md, encoding="utf-8")

    m6_files = sorted(
        p for p in out.rglob("*") if p.is_file() and p.name not in ("output_sha256.txt", ".gitkeep")
    )
    (out / "output_sha256.txt").write_text(
        "".join(f"{_sha256(p)}  {p.relative_to(ROOT)}\n" for p in m6_files),
        encoding="utf-8",
    )

    print(
        adjusted[
            [
                "cohort",
                "contrast",
                "lineage",
                "program",
                "n_patients",
                "mean_difference",
                "direction_status",
                "raw_p",
                "holm_adjusted_p",
            ]
        ].to_string(index=False)
    )


def self_check() -> None:
    """In-memory verification of N-candidate family Holm adjustment & 3-version comparison builder."""
    plan = pd.DataFrame(
        [
            {
                "test_id": f"{fam}_C{i}",
                "family": fam,
                "cohort": coh,
                "role": "external",
                "lineage": lin,
                "program": prog,
                "contrast": con,
                "expected_direction": "positive",
                "direction_rule": "mean_and_sign_majority",
                "alpha": 0.05,
                "test_method": "exact_sign_two_sided",
                "adjustment": "Holm",
                "analysis_tier": tier,
                "freeze_status": "frozen",
                "freeze_timestamp_utc": "2026-11-21T15:25:00+00:00",
                "freeze_version": "CRLM_V1.1_245_F1",
                "synthetic_example": False,
            }
            for i, (lin, prog) in enumerate(
                [("Mono/Mac/DC", "Core_Lipid_SAM_R"), ("T/NK", "CD8_Cytotoxic_Effector"), ("T/NK", "CD4_Treg_Suppressive")],
                start=1,
            )
            for coh, con, tier, fam in [
                ("GSE164522", "D_T", "primary", "PRIMARY_DT_164_V1_1"),
                ("GSE164522", "I", "background", "BACKGROUND_I_164_V1_1"),
                ("GSE178318", "D_T", "secondary", "SECONDARY_DT_178_V1_1"),
                ("GSE203608", "D_T", "secondary", "SECONDARY_DT_203608_V1_1"),
            ]
        ]
    )
    summaries = pd.DataFrame(
        [
            {
                "dataset": "GSE164522",
                "role": "external",
                "lineage": "Mono/Mac/DC",
                "program": "Core_Lipid_SAM_R",
                "contrast": "D_T",
                "n_patients": 8,
                "status": "estimated",
                "mean_difference": 0.045,
                "n_positive": 8,
                "n_negative": 0,
                "n_ties": 0,
                "n_sign_test": 8,
                "sign_p_two_sided": 2.0 / 256.0,
            }
        ]
    )
    adj, exp, fam, _ = adjust_registered(plan, summaries)
    assert len(adj) == 12 and fam["n_planned_tests"].eq(3).all()
    v1 = pd.read_csv(ARCHIVE_V1_M6) if ARCHIVE_V1_M6.is_file() else pd.DataFrame()
    tf = pd.read_csv(TECHFIX_M6) if TECHFIX_M6.is_file() else pd.DataFrame()
    comp_df, comp_md = build_three_version_comparison(v1, tf, adj, exp)
    assert not comp_df.empty and "CRLM_V1.1_245_F1" in comp_md
    print(f"assemble_validation --self-check PASSED: {len(comp_df)} comparison rows built (zero disk writes).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true", help="Run pure in-memory verification.")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        main()
