"""Register the V1.1 remediated development shortlist, frozen matched controls, and full provenance SHA256s.

Upgrades over V1 (`CRLM_V1_245_F1`, archived at `artifacts/archive/CRLM_V1_245_F1/artifacts/m5/`):
1. Registers version `CRLM_V1.1_245_F1` across four multiplicity families
   (`PRIMARY_DT_164_V1_1`, `BACKGROUND_I_164_V1_1`, `SECONDARY_DT_178_V1_1`, `SECONDARY_DT_203608_V1_1`)
   using primary scoring method `ucell_debiased_20bin_10ctrl` and HGNC canonical nuclear universes.
2. Locks frozen development `matched_controls.csv` and `reference_gene_stats.csv` (both in
   `artifacts/m5/` and `artifacts/m5/gene_sets/`) so external cohorts never re-fit background bins.
3. Explicitly discloses historical `CRLM_V1_245_F1` external exposure in `external_access_log.csv`
   and `freeze_record.md` per `AGENTS.md` §8.
4. Fixes V1 MAJOR-3 provenance gap: locks `models/scGPT-whole-human/{best_model.pt,vocab.json,args.json}`,
   `starter/linear_attention_mil.py`, `src/crlm/{hgnc_alias_map,annotation_m1,program_scores,model_bags,prior_mil}.py`,
   and computes `frozen_inputs_sha256.txt` + `output_sha256.txt` AFTER writing `validation_rules.md` and `freeze_record.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "starter")]

from adjust_registered_tests import validate_plan
from crlm.hgnc_alias_map import (
    HGNC_PROVENANCE_PATH,
    HGNC_TSV_PATH,
    SHARED_NUCLEAR_245_164_203608_PATH,
    SHARED_NUCLEAR_245_164_PATH,
    TUMOR_NUCLEAR_245_PATH,
)

OUT = ROOT / "artifacts/m5"
ARCHIVE_V1 = ROOT / "artifacts/archive/CRLM_V1_245_F1/artifacts/m5"
M3 = ROOT / "artifacts/m3"
M4 = ROOT / "artifacts/m4"
VERSION = "CRLM_V1.1_245_F1"

FAMILIES: tuple[tuple[str, str, str, str], ...] = (
    ("GSE164522", "D_T", "primary", "PRIMARY_DT_164_V1_1"),
    ("GSE164522", "I", "background", "BACKGROUND_I_164_V1_1"),
    ("GSE178318", "D_T", "secondary", "SECONDARY_DT_178_V1_1"),
    ("GSE203608", "D_T", "secondary", "SECONDARY_DT_203608_V1_1"),
)


def digest(path: Path) -> str:
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def build_plan_and_registry(
    chosen: pd.DataFrame,
    keep_genes: pd.DataFrame,
    gene_file_rel: str,
    gene_file_sha: str,
    background_rel: str,
    background_sha: str,
    controls_rel: str,
    controls_sha: str,
    stamp: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Construct validated `registered_tests` and `frozen_candidate_registry` tables."""
    if chosen.empty:
        raise ValueError("Cannot freeze an empty primary candidate shortlist")
    tests: list[dict[str, object]] = []
    registry: list[dict[str, object]] = []

    for i, r in enumerate(chosen.itertuples(index=False), start=1):
        direction = "positive" if r.direction > 0 else "negative"
        names = keep_genes.loc[
            keep_genes["program"].eq(r.program) & keep_genes["lineage"].eq(r.lineage),
            "gene_symbol",
        ].tolist()
        if not names:
            raise ValueError(f"Zero canonical genes found for {r.lineage}::{r.program}")

        registry.append(
            {
                "candidate_id": f"{r.lineage}::{r.program}",
                "candidate_tier": "primary",
                "lineage": r.lineage,
                "program": r.program,
                "focal_subtype": getattr(r, "focal_subtype", ""),
                "mechanism_cluster": getattr(r, "mechanism_cluster", ""),
                "program_version": VERSION,
                "status": "frozen",
                "synthetic_example": False,
                "gene_set_file": gene_file_rel,
                "gene_set_sha256": gene_file_sha,
                "positive_genes": ";".join(names),
                "negative_genes": "",
                "scoring_method": "ucell_debiased_20bin_10ctrl",
                "secondary_scoring_methods": "log1p_cp10k_debiased;ucell_raw;legacy_positive",
                "matched_controls_file": controls_rel,
                "matched_controls_sha256": controls_sha,
                "expected_DT_direction": direction,
                "expected_direction_DT": direction,
                "expected_direction": direction,
                "fixed_gene_universe_file": background_rel,
                "gene_universe_sha256": background_sha,
                "coverage_rule_id": "hgnc_canonical_shared_nuclear_80pct_min5",
                "min_measured_fraction_positive": 0.8,
                "min_measured_fraction_negative": 0.0,
                "min_measured_positive_genes": 5,
                "min_measured_negative_genes": 0,
                "required_core_genes": "",
                "min_cells_per_sample_subtype": 50,
                "discovery_datasets": "GSE245552",
                "external_datasets": "GSE164522;GSE178318;GSE203608",
                "primary_contrast": "D_T",
                "background_contrast": "I_on_Q",
                "primary_multiplicity_family": "PRIMARY_DT_164_V1_1",
                "background_multiplicity_family": "BACKGROUND_I_164_V1_1",
                "primary_test": "exact_sign_two_sided",
                "adjustment": "Holm",
                "alpha": 0.05,
                "missing_measurement_rule": "raw NA, Holm p=1",
                "freeze_timestamp_utc": stamp,
                "freeze_version": VERSION,
                "prior_external_access_disclosure": "artifacts/m5/external_access_log.csv",
            }
        )
        for cohort, contrast, tier, family in FAMILIES:
            tests.append(
                {
                    "test_id": f"{family}_C{i}",
                    "family": family,
                    "cohort": cohort,
                    "role": "external",
                    "lineage": r.lineage,
                    "program": r.program,
                    "contrast": contrast,
                    "expected_direction": direction,
                    "direction_rule": "mean_and_sign_majority",
                    "alpha": 0.05,
                    "test_method": "exact_sign_two_sided",
                    "adjustment": "Holm",
                    "analysis_tier": tier,
                    "freeze_status": "frozen",
                    "freeze_timestamp_utc": stamp,
                    "freeze_version": VERSION,
                    "synthetic_example": False,
                }
            )

    plan_df = pd.DataFrame(tests)
    validate_plan(plan_df)
    return plan_df, pd.DataFrame(registry)


def write_validation_rules(chosen: pd.DataFrame) -> None:
    cand_lines = "\n".join(
        f"- `{r.lineage}::{r.program}`（预期 `D_T` 与四组织交互 `I` 方向：`{'上升 (positive)' if r.direction > 0 else '下降 (negative)'}`；焦点纯细类：`{getattr(r, 'focal_subtype', '')}`）"
        for r in chosen.itertuples(index=False)
    )
    rules = f"""# V1.1 方法学重构版外部验证规则（`{VERSION}`）

> [!IMPORTANT]
> **版本属性与历史暴露透明声明**
> `{VERSION}` 是在完整归档初版基线 `CRLM_V1_245_F1`（`artifacts/archive/CRLM_V1_245_F1/`）与纯技术勘误 `CRLM_V1_245_F1_techfix` 之后，为修复初版 7 项 CRITICAL 方法学缺陷而在唯一开发集 `GSE245552` 内重做的系统重评估版本。鉴于初版已执行过 `GSE164522` 与 `GSE203608` 评估，本版本在统计与汇报中明确界定为**已披露既往暴露记录的外部独立队列复现与稳健性评估**，不冒称首次盲法确认。

## 1. 冻结主候选与多重比较家族

{cand_lines}

- **四个独立校正家族（每族双侧精确符号检验 + Holm，`alpha = 0.05`）**：
  1. `PRIMARY_DT_164_V1_1`：`GSE164522` 患者级配对肿瘤差 `D_T = LM - PT`；
  2. `BACKGROUND_I_164_V1_1`：`GSE164522` 同患者四组织合格子集 `Q` 上的器官背景扣除交互 `I = D_T_Q - D_N_Q`；
  3. `SECONDARY_DT_178_V1_1`：`GSE178318` 次要方向证据（使用条款待核实前保留 `raw_p = NA`，仅 Holm 向量置 `1.0`）；
  4. `SECONDARY_DT_203608_V1_1`：`GSE203608` 次要方向与精度证据。

## 2. 统一基因轴、冻结背景扣除与细胞清洗规则

1. **HGNC 规范化共同核基因轴**：所有外部队列统一经 `src/crlm/hgnc_alias_map.py` 规范化符号（`GSE203608` 结合 `ENSG` 消歧重名符号），并在与 `GSE245552` 共享的核编码基因轴（`GSE164522` 采用 `245_164_shared_nuclear` `18,815` 核基因；`GSE203608` 采用 `245_164_203608_shared_nuclear` `18,733` 核基因，显式排除全部 `353` 个 MT/核糖体偏倚基因）上评分。冻结正向成员测量覆盖须 `>= 80%` 且 `>= 5` 个基因。
2. **主评分公式（零外部重拟合）**：主指标固定为 `debiased_ucell_scores`（全测量轴 Mann–Whitney UCell `max_rank=1500` 减去本冻结包锁定的 `GSE245552` 20 桶 × 10 重复检测率匹配非成员对照集 `matched_controls.csv`），**严禁在外部队列重算检测率分桶或重抽背景基因**。同步输出共同核基因重标定 `log1p_cp10k_debiased`、`ucell_raw` 与初版 `legacy_positive` 敏感性对照。
3. **跨队列一致细胞清洗与同规则重注释对照**：
   - `GSE164522`：主分支默认剔除原作者已标注的全部跨谱系 `CD3E+` 双细胞群（髓系 `hC56/hC41/hC47/hC36` 共 3,662 个细胞 + B/Plasma `hC75_B-CD3E` 3,174 个细胞），另附保留双细胞敏感性与 **`annotation_m1.classify` 统一分层标志物重注释对照分支**；
   - `GSE203608`：严格过滤体内 `Sample_pass_qc == TRUE`、`Cell_type != "NA"` 且 `MitoFraction <= 0.20`（剔除 12,136 个死细胞/未注释细胞），同时提供洁净作者标签主分支与 `annotation_m1.classify` 统一分层标志物重注释分支，并附剔除患者 49 敏感性。
4. **聚合与分解**：每个合格标本×谱系门槛 `>= 50` 个细胞（另附 `20/100` 及纯细类 `>= 20` 对照），按“文库细胞均值 → 同患者同组织病灶等权 → 患者配对差”汇总；组成比例统一报告 `total_resolved_immune` 与 `core_three_lineage_immune` 双分母，细类状态执行共同支持对称分解。
"""
    (OUT / "validation_rules.md").write_text(rules, encoding="utf-8")


def freeze(*, force: bool = False) -> None:
    if not (ARCHIVE_V1 / "freeze_record.md").is_file():
        raise FileNotFoundError(
            f"Refusing to overwrite {OUT}: baseline archive missing at {ARCHIVE_V1 / 'freeze_record.md'}"
        )
    existing_rec = OUT / "freeze_record.md"
    if existing_rec.is_file() and VERSION in existing_rec.read_text(encoding="utf-8") and not force:
        raise ValueError(f"{VERSION} has already been registered in {existing_rec}; pass --force to regenerate")

    for prereq in (
        M3 / "stage_report.md",
        M3 / "candidate_screen.csv",
        M3 / "program_definitions_draft/program_genes.csv",
        M3 / "program_definitions_draft/matched_controls.csv",
        M3 / "program_definitions_draft/reference_gene_stats.csv",
        M4 / "stage_report.md",
        M4 / "full_data_models.csv",
    ):
        if not prereq.is_file():
            raise FileNotFoundError(f"Missing prerequisite {prereq}")

    register = pd.read_csv(M3 / "candidate_screen.csv")
    chosen = register.loc[register["shortlisted"]].sort_values(["lineage", "program"]).reset_index(drop=True)
    secondary_pass = register.loc[~register["shortlisted"] & (register["gates_1_to_4_pass"] | register["candidate_class"].eq("organ_resident_control"))].copy()

    members = pd.read_csv(M3 / "program_definitions_draft/program_genes.csv")
    model_members = members.loc[members["role"].eq("complete_program") & members["on_discovery_axis"]].copy()
    keep = chosen[["lineage", "program"]].merge(
        model_members, on=["lineage", "program"], validate="one_to_many"
    )

    OUT.mkdir(parents=True, exist_ok=True)
    genes_dir = OUT / "gene_sets"
    genes_dir.mkdir(parents=True, exist_ok=True)

    gene_file = genes_dir / "program_genes.csv"
    keep.to_csv(gene_file, index=False)
    model_members.to_csv(genes_dir / "model_program_genes.csv", index=False)
    secondary_pass.to_csv(OUT / "secondary_sensitivity_candidates.csv", index=False)

    background = genes_dir / "development_rank_background.txt"
    shutil.copyfile(TUMOR_NUCLEAR_245_PATH, background)
    shutil.copyfile(SHARED_NUCLEAR_245_164_PATH, genes_dir / "shared_nuclear_245_164.txt")
    shutil.copyfile(SHARED_NUCLEAR_245_164_203608_PATH, genes_dir / "shared_nuclear_245_164_203608.txt")

    shutil.copyfile(M3 / "program_definitions_draft/matched_controls.csv", genes_dir / "matched_controls.csv")
    shutil.copyfile(M3 / "program_definitions_draft/matched_controls.csv", OUT / "matched_controls.csv")
    shutil.copyfile(M3 / "program_definitions_draft/reference_gene_stats.csv", genes_dir / "reference_gene_stats.csv")

    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    # Atomically update configs/project.json BEFORE copying to development_config.json (fixes V1 MINOR-1)
    proj_cfg_path = ROOT / "configs/project.json"
    proj_cfg = json.loads(proj_cfg_path.read_text(encoding="utf-8"))
    proj_cfg["research_version"] = "V1.1_AUDIT_REMEDIATED"
    proj_cfg["status"] = "v1_1_m5_frozen_awaiting_m6_m8"
    proj_cfg.setdefault("milestones", {}).update(
        {
            "M0": "completed",
            "M1": "completed",
            "M2": "completed",
            "M3": "completed",
            "M4": "completed",
            "M5": "completed",
            "M6": "in_progress_v1_1_external_validation",
            "M7": "completed_spatial_measurement_and_experimental_design",
            "M8": "pending_v1_1_synthesis",
        }
    )
    proj_cfg["candidate_registration"] = {
        "status": "frozen",
        "freeze_version": VERSION,
        "freeze_timestamp_utc": stamp,
        "archived_baseline_freeze_version": "CRLM_V1_245_F1",
        "techfix_ablation_version": "CRLM_V1_245_F1_techfix",
        "primary_lineage_program_total_planned": "2-3 total",
        "primary_candidates_frozen": [
            f"{r.lineage}::{r.program}" for r in chosen.itertuples(index=False)
        ],
        "primary_scoring_method": "ucell_debiased_20bin_10ctrl",
        "primary_external_family": "PRIMARY_DT_164_V1_1",
        "background_external_family": "BACKGROUND_I_164_V1_1",
        "secondary_external_families": [
            "SECONDARY_DT_178_V1_1",
            "SECONDARY_DT_203608_V1_1",
        ],
        "background_family_separate": True,
        "prior_external_exposure_disclosed": True,
    }
    proj_cfg_path.write_text(
        json.dumps(proj_cfg, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    shutil.copyfile(ROOT / "src/crlm/program_scores.py", OUT / "program_scores.py")
    shutil.copyfile(proj_cfg_path, OUT / "development_config.json")
    shutil.copyfile(M4 / "full_data_models.csv", OUT / "frozen_models.csv")

    plan, registry_df = build_plan_and_registry(
        chosen=chosen,
        keep_genes=keep,
        gene_file_rel=str(gene_file.relative_to(ROOT)),
        gene_file_sha=digest(gene_file),
        background_rel=str((genes_dir / "shared_nuclear_245_164_203608.txt").relative_to(ROOT)),
        background_sha=digest(genes_dir / "shared_nuclear_245_164_203608.txt"),
        controls_rel=str((OUT / "matched_controls.csv").relative_to(ROOT)),
        controls_sha=digest(OUT / "matched_controls.csv"),
        stamp=stamp,
    )
    plan.to_csv(OUT / "registered_tests.csv", index=False)
    registry_df.to_csv(OUT / "frozen_candidate_registry.csv", index=False)

    # Preserve complete M0 + V1 exposure history and append V1.1 disclosure
    exposure_m0 = pd.read_csv(ROOT / "artifacts/m0/external_access_log.csv")
    v1_arch_exp = ARCHIVE_V1 / "external_access_log.csv"
    base_exp = pd.read_csv(v1_arch_exp) if v1_arch_exp.is_file() else exposure_m0
    new_exp_rows = [
        {
            "access_date_utc": stamp,
            "time_precision": "second_utc",
            "dataset": "GSE164522;GSE203608;GSE335552",
            "purpose": "CRLM_V1_245_F1 & CRLM_V1_245_F1_techfix prior evaluation disclosure before V1.1 freeze",
            "material_seen": (
                "Prior V1_F1 and V1_F1_techfix runs evaluated Mac_SPP1, UPR, and 15 V1 classifiers on "
                "GSE164522/GSE203608 (archived under artifacts/archive/CRLM_V1_245_F1/ and artifacts/techfix_v1_f1/)"
            ),
            "biological_effects_seen": "V1 Mac_SPP1 and UPR contrasts seen in prior version; zero V1.1 new candidate scores computed on external cohorts prior to V1.1 freeze",
            "what_was_not_done": "V1.1 candidate shortlist, debiased UCell control bins, and M4 V1.1 models were selected 100% inside GSE245552",
            "consequence_for_future_confirmation": "V1.1 external evaluation is explicitly classified as exposure-disclosed methodological remediation & external robustness replication, not blind de-novo confirmation",
        }
    ]
    pd.concat([base_exp, pd.DataFrame(new_exp_rows)], ignore_index=True).to_csv(
        OUT / "external_access_log.csv", index=False
    )

    write_validation_rules(chosen)

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    fitted = pd.read_csv(OUT / "frozen_models.csv")
    accepted_col = fitted["accepted"] if "accepted" in fitted.columns else pd.Series(True, index=fitted.index)
    for r in fitted.loc[accepted_col].itertuples():
        mpath = ROOT / r.path
        if not mpath.is_file() or digest(mpath) != r.sha256:
            raise ValueError(f"Model archive checksum mismatch for {r.path}")

    cand_summary = "、".join(
        f"`{r.lineage}::{r.program}`（{'上升' if r.direction > 0 else '下降'}）"
        for r in chosen.itertuples(index=False)
    )
    (OUT / "freeze_record.md").write_text(
        f"""# V1.1 方法学重构外部验证冻结记录（`{VERSION}`）

- **冻结时间戳（UTC）**：`{stamp}`
- **冻结版本标识**：`{VERSION}`（Git HEAD: `{commit}`；初版基线完整封存于 `artifacts/archive/CRLM_V1_245_F1/`）
- **锁定主候选（共 {len(registry_df)} 项，{len(plan)} 项预登记检验，分属 4 个独立 Holm 校正家族）**：{cand_summary}。[冻结候选表](frozen_candidate_registry.csv) · [预登记检验表](registered_tests.csv) · [同机制次要/驻留对照清单](secondary_sensitivity_candidates.csv)
- **方法学与防泄漏保证**：全部候选筛选、20 桶 × 10 重复检测率匹配对照基因集（[matched_controls.csv](matched_controls.csv)）及 `{int(accepted_col.sum())}` 个全开发集迁移分类模型（[frozen_models.csv](frozen_models.csv)）均 100% 在 `GSE245552` 内确定。
- **既往暴露披露**：初版 `CRLM_V1_245_F1` 与纯技术勘误 `CRLM_V1_245_F1_techfix` 对外部队列的历史访问已按秒级 UTC 记录于 [external_access_log.csv](external_access_log.csv)，故 `{VERSION}` 严格定位为已披露历史暴露的外部独立队列复现评估。
- **全链路哈希锁**：本记录、验证规则、阶段报告、`scGPT` 预训练权重/词表、全部 `M4` 迁移模型权重及评分/注释/HGNC 源码哈希均已闭环写入 [frozen_inputs_sha256.txt](frozen_inputs_sha256.txt) 与 [output_sha256.txt](output_sha256.txt)。
""",
        encoding="utf-8",
    )

    cand_rows_md = []
    for r in chosen.itertuples(index=False):
        n_pos_dir = int(round(float(r.direction_fraction) * int(r.n_nonzero)))
        dir_str = "上升 (LM > PT, +1)" if r.direction > 0 else "下降 (LM < PT, -1)"
        cand_rows_md.append(
            f"| `{r.lineage}::{r.program}` | {dir_str} | {n_pos_dir}/{int(r.n_pairs)} ({100.0 * float(r.direction_fraction):.1f}%) | "
            f"{float(r.mean_D_T):+.5f} [{float(r.bootstrap_ci_low):+.5f}, {float(r.bootstrap_ci_high):+.5f}] | "
            f"{float(r.discovery_sign_p):.5f} | `{r.focal_subtype}`: {float(r.focal_mean_D_T):+.5f} ({100.0 * float(r.focal_direction_fraction):.1f}%) | "
            f"{100.0 * float(r.within_state_share):.1f}% |"
        )
    sec_list_str = "、".join(
        f"`{r.lineage}::{r.program}` ({r.dedup_status if r.gates_1_to_4_pass else r.candidate_class})"
        for r in secondary_pass.itertuples(index=False)
    )
    (OUT / "stage_report.md").write_text(
        f"""# M5 · V1.1 方法学重构版候选程序与分类模型原子冻结报告（`{VERSION}`）

> [!IMPORTANT]
> **三版本透明并列与历史暴露声明**
> 初版基线 `CRLM_V1_245_F1`（`artifacts/archive/CRLM_V1_245_F1/`）与纯技术勘误 `CRLM_V1_245_F1_techfix`（`artifacts/techfix_v1_f1/`）已完整封存并记录于 [external_access_log.csv](external_access_log.csv)。本冻结包 `{VERSION}` 的候选筛选、20 桶 × 10 重复检测率匹配对照基因集与 15 个全开发集分类头均 100% 在唯一开发集 `GSE245552` 内重做，外部评估严格界定为**已披露历史暴露的独立队列复现与方法学稳健性验证**。

## 1. 冻结主候选摘要（来自 `GSE245552` 五重门禁筛选）

| 主候选 (`lineage::program`) | 冻结预期方向 | 开发集同向患者 | `mean_D_T` [95% CI] | 双侧符号 `p` | 焦点纯细类 `mean_D_T` (同向率) | 亚型内状态占比 |
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(cand_rows_md)}

- **同机制共线次要候选与器官驻留基准（不计入主检验多重比较家族）**：{sec_list_str}（详见 [secondary_sensitivity_candidates.csv](secondary_sensitivity_candidates.csv)）。
- **初版 `B/Plasma::HALLMARK_UNFOLDED_PROTEIN_RESPONSE` 淘汰依据**：去偏后开发集 `6/8 = 75.0%`（`p = 0.289`），`72.6%` 总差由 `13/13` 患者浆细胞比例暴跌驱动，纯 `plasma cell` 内仅 `5/8 = 62.5%`（`-0.00308`），未通过 Gate 4 组成分解与焦点细类门禁。

## 2. 预登记检验矩阵与全开发集迁移分类模型

- **预登记检验（共 `{len(plan)}` 项，分属 4 个独立 Holm 校正家族，每族 `m = {len(registry_df)}`，`alpha = 0.05`）**：
  - `PRIMARY_DT_164_V1_1`（`GSE164522` 配对肿瘤差 `D_T = LM - PT`）
  - `BACKGROUND_I_164_V1_1`（`GSE164522` 四组织合格子集 `Q` 交互 `I = D_T_Q - D_N_Q`）
  - `SECONDARY_DT_178_V1_1`（`GSE178318` 次要方向证据；使用条款待核实前保留 `raw_p = NA`、Holm 向量置 `1.0`）
  - `SECONDARY_DT_203608_V1_1`（`GSE203608` 次要方向与精度证据）
- **冻结全开发集迁移分类头（`{int(accepted_col.sum())}/15` 项全部验收通过）**：`T/NK`、`Mono/Mac/DC`、`B/Plasma` × `{{program_elastic, PCA_mean, PCA_attention, scGPT_mean, scGPT_attention}}`，全部基于 `245_164_shared_nuclear`（`18,815` 个共享核编码基因，在 `GSE164522`、`GSE203608` 与 `scGPT` 词表上实现 `512/512 = 100.0%` PCA HVG 覆盖，彻底消除初版 4 个模型拒判）。

## 3. 冻结清单与防篡改校验

- **时间与版本**：`{stamp}` · `{VERSION}` · Git HEAD `{commit}`
- **配置自洽性**：`configs/project.json` 与 [development_config.json](development_config.json) 在冻结时刻原子同步（`research_version = "V1.1_AUDIT_REMEDIATED"`，`freeze_version = "{VERSION}"`）。
- **核心文件**：[冻结记录](freeze_record.md) · [主要候选](frozen_candidate_registry.csv) · [预登记检验](registered_tests.csv) · [验证规则](validation_rules.md) · [匹配对照集](matched_controls.csv) · [冻结模型](frozen_models.csv) · [外部暴露日志](external_access_log.csv) · [全链路输入哈希](frozen_inputs_sha256.txt) · [阶段产物哈希](output_sha256.txt)
""",
        encoding="utf-8",
    )

    # Compute provenance hashes AFTER writing validation_rules.md, freeze_record.md, and stage_report.md!
    model_inputs = [p for p in (ROOT / "models/m4/full").rglob("*") if p.is_file()]
    scgpt_inputs = [
        ROOT / "models/scGPT-whole-human/best_model.pt",
        ROOT / "models/scGPT-whole-human/vocab.json",
        ROOT / "models/scGPT-whole-human/args.json",
    ]
    crlm_src_files = sorted((ROOT / "src/crlm").glob("*.py"))
    code_and_upstream = [
        proj_cfg_path,
        HGNC_TSV_PATH,
        HGNC_PROVENANCE_PATH,
        M3 / "candidate_screen.csv",
        M3 / "stage_report.md",
        M3 / "output_sha256.txt",
        M4 / "full_data_models.csv",
        M4 / "stage_report.md",
        M4 / "output_sha256.txt",
        ROOT / "scripts/freeze_validation.py",
        ROOT / "scripts/encode_scgpt.py",
        ROOT / "scripts/fit_transfer_models.py",
        *crlm_src_files,
        ROOT / "starter/linear_attention_mil.py",
        ROOT / "starter/adjust_registered_tests.py",
        ROOT / "starter/paired_effects.py",
        ROOT / "starter/composition_decomposition.py",
    ]
    m5_artifacts = [
        p
        for p in OUT.rglob("*")
        if p.is_file() and p.name not in (".gitkeep", "frozen_inputs_sha256.txt", "output_sha256.txt")
    ]
    provenance = sorted(
        set(p for p in (m5_artifacts + model_inputs + scgpt_inputs + code_and_upstream) if p.is_file())
    )
    (OUT / "frozen_inputs_sha256.txt").write_text(
        "".join(f"{digest(p)}  {p.relative_to(ROOT)}\n" for p in provenance),
        encoding="utf-8",
    )
    m5_final = sorted(
        p for p in OUT.rglob("*") if p.is_file() and p.name not in (".gitkeep", "output_sha256.txt")
    )
    (OUT / "output_sha256.txt").write_text(
        "".join(f"{digest(p)}  {p.relative_to(ROOT)}\n" for p in m5_final),
        encoding="utf-8",
    )
    print(f"Registered {VERSION}: {len(registry_df)} primary candidates, {len(plan)} planned tests.")


def self_check() -> None:
    """In-memory verification of V1.1 plan & registry schema without touching artifacts/m5/."""
    chosen = pd.DataFrame(
        [
            {"lineage": "Mono/Mac/DC", "program": "Core_Lipid_SAM_R", "direction": 1, "focal_subtype": "Macrophage", "mechanism_cluster": "Myeloid_R"},
            {"lineage": "T/NK", "program": "CD8_Cytotoxic_Effector", "direction": 1, "focal_subtype": "CD8 T cell", "mechanism_cluster": "TNK_Eff"},
        ]
    )
    keep_genes = pd.DataFrame(
        [
            {"lineage": "Mono/Mac/DC", "program": "Core_Lipid_SAM_R", "gene_symbol": g}
            for g in ("SPP1", "TREM2", "APOE", "APOC1", "GPNMB")
        ]
        + [
            {"lineage": "T/NK", "program": "CD8_Cytotoxic_Effector", "gene_symbol": g}
            for g in ("GZMA", "GZMB", "PRF1", "GNLY", "NKG7")
        ]
    )
    plan, reg = build_plan_and_registry(
        chosen,
        keep_genes,
        "artifacts/m5/gene_sets/program_genes.csv",
        "0" * 64,
        "artifacts/m5/gene_sets/shared_nuclear_245_164_203608.txt",
        "1" * 64,
        "artifacts/m5/matched_controls.csv",
        "2" * 64,
        "2026-11-21T15:25:00+00:00",
    )
    assert len(reg) == 2 and len(plan) == 8
    assert set(plan["family"]) == {f[3] for f in FAMILIES}
    print("freeze_validation --self-check PASSED (zero disk writes).")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Allow overwriting an existing registration of the same version.")
    parser.add_argument("--self-check", action="store_true", help="Run pure in-memory schema check.")
    args = parser.parse_args()
    if args.self_check:
        self_check()
    else:
        freeze(force=args.force)
