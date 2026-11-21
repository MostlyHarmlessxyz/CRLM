# 复算入口

所有命令从项目根目录运行。各阶段锁定的 Python 依赖见阶段报告。原始表达、模型权重、H5AD 及表示缓存保存在各自数据目录，Git 保存来源、变换、输入标识、分析表和图。已有发布版本复核时使用 Git 中的冻结包 `72ee86b`。

| 步骤 | 入口及依赖 | 复算输出 |
|---|---|---|
| 环境与身份 | [M0 原命令与隔离自检](../m0/stage_report.md) | 真实患者映射、分折、历史与本机检查 |
| 开发集获取和导入 | `scripts/fetch_gse245552.py`、[Pt18 导入命令](../m1/stage_report.md)、[剩余文库与统一处理](../m2/stage_report.md) | 39 个对象、QC、谱系覆盖和患者伪批量 |
| 程序与计数 | [M3 三组脚本与环境](../m3/stage_report.md) | 患者差、分解、技术敏感性和计数模型 |
| 分类 | [M4 顺序与 A750 实测](../m4/stage_report.md)、[权重准备](../m4/resource_profiles/scgpt_sources.md)、[直接前作](../m4/direct_prior/stage_report.md) | 折外预测、归因、15 个迁移头 |
| 冻结 | `sha256sum -c artifacts/m5/frozen_inputs_sha256.txt` | 已登记的判据、表达成员、程序和模型 |
| 主要外部来源 | `scripts/validate_gse164522.py fetch`，随后依次运行 `identity`、`convert --freeze-commit 72ee86b`、`score --freeze-commit 72ee86b` | 四组织处理表达对象、标本与患者程序效应 |
| 次要外部来源 | `scripts/fetch_gse203608.py` 和 `scripts/validate_gse203608.py` | 原版本体内计数核验、四位解剖配对患者及谱系有效 N |
| 冻结迁移预测 | [输入尺度](../m6/GSE164522/model_predictions/normalization_source.md)，按预测入口的 `development-check`、`encode`、`predict` 运行，外部调用给定 `--freeze-commit 72ee86b` | 相同肿瘤细胞下的固定头预测、覆盖、区间 |
| 净效应与校正 | `scripts/assemble_validation.py` | 同一患者的四组织交互核对、全部八项登记检验 |
| 局部空间 | [空间输入和规则](../m7/spatial_measurement_audit/spatial_rules.md)，运行 `scripts/spatial_support.py` | 患者、切片及实测标记距离 |
| 主图 | `scripts/build_paper_figures.py` | PNG、矢量 PDF 及其源数据 |

GSE335552 可按 GEO 的原始附件表取得所需八个文件。以下命令只选细胞矩阵和质心表，在冷启动工作区中重建本研究实际使用的原件。

```bash
mkdir -p data/raw/GSE335552
awk '/^!Sample_supplementary_file = / && /(_cell_feature_matrix.h5|_cells.parquet.gz)$/ {sub("^ftp:", "https:", $3); print $3}' starter/audit_data/geo/GSE335552_metadata_only.soft |
while IFS= read -r url; do
  curl -fLsS --retry 3 -o "data/raw/GSE335552/${url##*/}" "$url" || exit 1
done
```
