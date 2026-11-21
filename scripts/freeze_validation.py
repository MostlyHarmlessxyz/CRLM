"""Register the development shortlist and its validation families before external scoring."""

import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "starter"))
from adjust_registered_tests import validate_plan

OUT = ROOT / "artifacts/m5"
M3 = ROOT / "artifacts/m3"
M4 = ROOT / "artifacts/m4"
VERSION = "CRLM_V1_245_F1"


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def freeze():
    if (OUT / "freeze_record.md").exists():
        raise ValueError("This validation version has already been registered")
    for report in (
        M3 / "stage_report.md",
        M4 / "stage_report.md",
        OUT / "validation_rules.md",
    ):
        if not report.is_file():
            raise ValueError(f"Missing prerequisite {report}")
    register = pd.read_csv(M3 / "candidate_screen.csv")
    chosen = register[register.shortlisted].sort_values(["lineage", "program"])
    members = pd.read_csv(M3 / "program_definitions_draft/program_genes.csv")
    model_members = members[
        members.role.eq("complete_program") & members.on_discovery_axis
    ].copy()
    keep = chosen[["lineage", "program"]].merge(
        model_members, on=["lineage", "program"], validate="one_to_many"
    )
    genes_dir = OUT / "gene_sets"
    genes_dir.mkdir(parents=True, exist_ok=True)
    gene_file = genes_dir / "program_genes.csv"
    keep.to_csv(gene_file, index=False)
    model_members.to_csv(genes_dir / "model_program_genes.csv", index=False)
    background = genes_dir / "development_rank_background.txt"
    shutil.copyfile(M3 / "program_definitions_draft/ranking_background.txt", background)
    code = OUT / "program_scores.py"
    shutil.copyfile(ROOT / "src/crlm/program_scores.py", code)
    shutil.copyfile(ROOT / "configs/project.json", OUT / "development_config.json")
    shutil.copyfile(M4 / "full_data_models.csv", OUT / "frozen_models.csv")
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    tests = []
    registry = []
    for i, r in enumerate(chosen.itertuples(index=False), start=1):
        direction = "positive" if r.direction > 0 else "negative"
        names = keep.loc[
            keep.program.eq(r.program) & keep.lineage.eq(r.lineage), "gene_symbol"
        ]
        registry.append(
            dict(
                candidate_id=r.lineage + "::" + r.program,
                lineage=r.lineage,
                program=r.program,
                program_version=VERSION,
                status="frozen",
                synthetic_example=False,
                gene_set_file=str(gene_file.relative_to(ROOT)),
                gene_set_sha256=digest(gene_file),
                positive_genes=";".join(names),
                negative_genes="",
                scoring_method="positive_midrank_1500",
                expected_DT_direction=direction,
                expected_direction=direction,
                fixed_gene_universe_file=str(background.relative_to(ROOT)),
                gene_universe_sha256=digest(background),
                coverage_rule_id="source_common_axis_exact_symbol_80pct_min5",
                min_measured_fraction_positive=0.8,
                min_measured_fraction_negative=0,
                min_measured_positive_genes=5,
                min_measured_negative_genes=0,
                required_core_genes="",
                min_cells_per_sample_subtype=50,
                discovery_datasets="GSE245552",
                external_datasets="GSE164522;GSE178318;GSE203608",
                primary_contrast="D_T",
                background_contrast="I_on_Q",
                primary_multiplicity_family="PRIMARY_DT_164_V1",
                background_multiplicity_family="BACKGROUND_I_164_V1",
                primary_test="exact_sign_two_sided",
                adjustment="Holm",
                alpha=0.05,
                missing_measurement_rule="raw NA, Holm p=1",
                freeze_timestamp_utc=stamp,
                freeze_version=VERSION,
                prior_external_access_disclosure="artifacts/m5/external_access_log.csv",
            )
        )
        for cohort, contrast, tier, family in [
            ("GSE164522", "D_T", "primary", "PRIMARY_DT_164_V1"),
            ("GSE164522", "I", "background", "BACKGROUND_I_164_V1"),
            ("GSE178318", "D_T", "secondary", "SECONDARY_DT_178_V1"),
            ("GSE203608", "D_T", "secondary", "SECONDARY_DT_203608_V1"),
        ]:
            tests.append(
                dict(
                    test_id=f"{family}_C{i}",
                    family=family,
                    cohort=cohort,
                    role="external",
                    lineage=r.lineage,
                    program=r.program,
                    contrast=contrast,
                    expected_direction=direction,
                    direction_rule="mean_and_sign_majority",
                    alpha=0.05,
                    test_method="exact_sign_two_sided",
                    adjustment="Holm",
                    analysis_tier=tier,
                    freeze_status="frozen",
                    freeze_timestamp_utc=stamp,
                    freeze_version=VERSION,
                    synthetic_example=False,
                )
            )
    plan = pd.DataFrame(tests)
    validate_plan(plan)
    plan.to_csv(OUT / "registered_tests.csv", index=False)
    pd.DataFrame(registry).to_csv(OUT / "frozen_candidate_registry.csv", index=False)
    exposure = pd.read_csv(ROOT / "artifacts/m0/external_access_log.csv")
    rows = [
        dict(
            access_date_utc=stamp[:10],
            time_precision="date_only",
            dataset=g,
            purpose="Pre-freeze technical preparation",
            material_seen=material,
            biological_effects_seen="No new program effects",
            what_was_not_done="No new candidate scores or contrasts",
            consequence_for_future_confirmation="See historical access rows",
        )
        for g, material in [
            (
                "GSE164522",
                "Author identity metadata, tissue and sample coding; compressed files staging",
            ),
            (
                "GSE203608",
                "Original identity metadata, gene annotation, MatrixMarket header and compressed staging",
            ),
            (
                "GSE335552",
                "Archived measured feature lists, H5 axes and centroid schema; public files staging",
            ),
            ("GSE178318", "Archived GEO usage terms; no new expression access"),
        ]
    ]
    pd.concat([exposure, pd.DataFrame(rows)]).to_csv(
        OUT / "external_access_log.csv", index=False
    )
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    fitted = pd.read_csv(OUT / "frozen_models.csv")
    if not all(
        (ROOT / r.path).is_file() and digest(ROOT / r.path) == r.sha256
        for r in fitted[fitted.accepted].itertuples()
    ):
        raise ValueError("Model archive does not match registered checksums")
    model_inputs = list((ROOT / "models/m4/full").glob("*/*"))
    provenance = list(
        sorted(
            set(
                [p for p in OUT.rglob("*") if p.is_file() and p.name != ".gitkeep"]
                + model_inputs
                + [
                    M3 / "candidate_screen.csv",
                    M4 / "stage_report.md",
                    ROOT / "scripts/freeze_validation.py",
                    ROOT / "scripts/encode_scgpt.py",
                    ROOT / "src/crlm/model_bags.py",
                    ROOT / "scripts/fit_transfer_models.py",
                    ROOT / "starter/adjust_registered_tests.py",
                    ROOT / "starter/paired_effects.py",
                ]
            )
        )
    )
    (OUT / "frozen_inputs_sha256.txt").write_text(
        "".join(f"{digest(p)}  {p.relative_to(ROOT)}\n" for p in provenance)
    )
    (OUT / "freeze_record.md").write_text(f"""# V1 外部验证冻结

实际冻结时间为 {stamp}，版本为 `{VERSION}`。开发源码与阶段结果的 Git 版本为 `{commit}`。本次将 {len(registry)} 项候选、{len(plan)} 项登记检验锁定为四个各自校正的家族。[候选表](frozen_candidate_registry.csv) · [登记检验](registered_tests.csv)

主候选为髓系 Mac_SPP1 上升和 B/Plasma UPR 下降。后者在发现集对深度与组成敏感。所有登记项均须保留，次要来源 178 的使用条款另行落实。此前 164 已发表亚群比例的暴露已登记。[验证规则](validation_rules.md) · [访问记录](external_access_log.csv)

单细胞秩函数、成员及发现排序轴已归档。各目标来源按事先指定规则生成其自身的共同测量轴。分类压力测试绑定 {int(fitted.accepted.sum())} 个全开发集拟合的模型及训练标本。[模型清单](frozen_models.csv) · [源码、输入与模型 SHA](frozen_inputs_sha256.txt)
""")


if __name__ == "__main__":
    freeze()
