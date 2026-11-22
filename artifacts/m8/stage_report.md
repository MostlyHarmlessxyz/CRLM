# M8 · `CRLM_V1.1_AUDIT_REMEDIATED` 终稿综合、出版级主图与三版本可复现交付报告

1. **三版本透明并列终稿与四幅主图**：
   - 已从最新阶段结果表正式运行 `scripts/build_paper_figures.py`，生成四幅高分辨率主图（PNG + 矢量 PDF）与 7 份逐图对应源数据表（`source_data/`）：
     - `01_patient_programs.{png,pdf}`：展示 `V1.1` 两个新主候选（`Mono/Mac/DC::Core_Lipid_SAM_R` 与 `T/NK::B07_Th17_TWEAK`）在 `GSE245552`、`GSE164522`、`GSE203608` 的逐患者配对差散点与均值符号计数、**`GSE178318`（`N=6` 未暴露 `terms_pending`）显式预登记占位面板**、`GSE164522` 同患者四组织子集 `Q`（`D_T_Q / D_N_Q / I` 逐供者连线与 `1+/9−` 分布），以及 `B/Plasma::UPR` 对称细类分解与四组织近似抵消证伪对照；
     - `02_development_classifiers.{png,pdf}`：并列呈现 `GSE245552` 单谱系全袋模型（含修复后的 `PCA_subtype_balanced`、折内实质基因掩码与 `200` 次患者置换 95% 经验零分布带）与 `All3Lineages` 跨三谱系联合标本模型；
     - `03_spatial_distances.{png,pdf}`：呈现 `GSE335552` 接触态 `CD8 T` 救回后的 `400 -> 200 -> 100 -> 50 µm` 多尺度距离差变化曲线与近 `EPCAM+` 肿瘤上皮 `<15 µm` vs 间质带 `>=15 µm` 分层对比（清晰标注各切片 `CD8 T` 数量与合格方格数）；
     - `04_external_classifiers.{png,pdf}`：并列呈现 `GSE164522` 外部迁移在 `V1.1`（`15/15` 全袋零拒判，含 95% 区间）、`V1_F1_techfix`（`15/15`，含 95% 区间）与归档基线 `V1_F1`（`11/15` 点估计标记及 4 处显式 `<95%` HVG 拒判标注）的三版本患者组织级 AUROC（对应 `V1.1` 源表见 `source_data/external_classifier_metrics.csv`，三版本并列表见 `../m6/three_version_classifier_comparison.csv`）。
   - [最终实证报告](final_report.md) · [局限与方案偏差](limitations_and_deviations.md) · [图目录](figures/) · [图源数据](source_data/)

2. **全链路复算清单与路径 100% 实证核验**：
   - [reproducibility_manifest.csv](reproducibility_manifest.csv) 完整登记 `M0–M8` 与 `Techfix` 共 24 条计算步骤，经回归测试脚本逐一核实全部 120 个文件与路径条目真实存在、全部 24 个隔离 `.venv-*` Python 解释器具备可执行权限（`+x`），彻底消除了初版路径与参数错填问题。[复算指南](reproduce.md)

3. **阶段哈希闭环校验**：
   ```bash
   sha256sum -c artifacts/m8/output_sha256.txt
   ```
   根级入口 `STATUS.md` 与 `README.md` 已同步至 `CRLM_V1.1_AUDIT_REMEDIATED` 实证状态，明确保留 `GSE178318`（待 PI 使用条件）、`GSE225857`（缺供者映射排除）及人源组织原位与同供者趋化验证预案边界。
