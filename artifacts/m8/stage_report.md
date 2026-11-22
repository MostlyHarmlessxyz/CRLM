# M8 · `CRLM_V1.1_AUDIT_REMEDIATED` 终稿综合、8/8 全观测外部闭环与双队列主图交付报告

1. **8/8 预登记检验全观测闭环与四幅确定性出版级主图（PNG + 矢量 PDF + 9 份逐元素全等源表）**：
   - 依托 M5（`5b01467`）零字节改动的冻结契约完成 `GSE178318` 前瞻零表达暴露验证与真实三版本子袋 XPU 编码后，已重跑 `scripts/build_paper_figures.py`（固定 PDF `CreationDate=None` 与 PNG 软件元数据，重跑字节级幂等）更新四幅主图与全部 **9 份逐元素全等源表**（`source_data/`，新增 `three_version_classifier_comparison.csv` 覆盖图 4 全部 `90` 条三版本分类记录）：
     - `01_patient_programs.{png,pdf}`：完整绘出 `Core_Lipid_SAM_R` 与 `B07_Th17_TWEAK` 在全部四个计划队列（`245` 开发、`164` 主外部、`203608` 次要及**淡绿底高亮的 `178 (Zero-Exp, N=6)` 零表达暴露前瞻队列**）的逐患者配对散点、`GSE164522` 同患者四组织 `Q` 子集（`B07_Th17_TWEAK` 观测交互 `I = −0.04918, Holm p = 0.0430 <= 0.05`）以及三队列（`245 / 164 / 178`，注明 `203608 N/A` 及 `178 N=4: 65.6% Cmp, 34.4% St`）`B/Plasma::UPR` 对称细类分解对照；
     - `02_development_classifiers.{png,pdf}`：单谱系全袋模型（含 `PCA_subtype_balanced = 0.990`、实质基因掩码与 `200-perm 95% null envelope [0.22, 0.77]` 外包络带）与 `All3Lineages` 跨三谱系联合标本模型；
     - `03_spatial_distances.{png,pdf}`：`GSE335552` 救回接触态 `CD8 T` 后的 `400 -> 200 -> 100 -> 50 µm` 多尺度衰减曲线与 `<15 µm EPCAM+` 肿瘤巢 vs `>=15 µm` 间质带对称分层对比；
     - `04_external_classifiers.{png,pdf}`：`2 × 3` 双外部队列面板，上排为 `GSE164522`（暴露后重评估），下排为 **`GSE178318`（零表达暴露前瞻验证）**，并列展示 `V1.1`（`15/15` 全袋零拒判）、`V1_F1_techfix`（`15/15`，独立 HGNC 分词子袋 XPU 编码）与原版 `V1_F1`（两队列重现一致的 `4/15` `<95%` 原始符号 HVG 拒判，独立原始轴子袋 XPU 编码），呈现零样本 `scGPT` 未相对简约线性 `PCA` 提供一致外部增量（`178` `T/NK` `PCA_mean − scGPT_mean = +0.333 [+0.083, +0.556]`，而 `10/12` 对 CI 跨 `0` 且 `B/Plasma` 注意力 `scGPT 0.667 > PCA 0.625`）。
   - [最终实证报告](final_report.md) · [局限与方案偏差](limitations_and_deviations.md) · [图目录](figures/) · [图源数据](source_data/)

2. **26 行全链路直接生产脚本复算索引与路径 100% 核实**：
   - [reproducibility_manifest.csv](reproducibility_manifest.csv) 严格按直接生产脚本、直接产物文件与收归仓内的可执行 `.venv-*` 解释器登记 `M0–M8`（含 `starter/validate_manifest.py`、`fetch_gse178318.py`、`validate_gse178318.py`、`predict_gse164522.py` 双阶段、`predict_gse178318.py`）与 `Techfix` 分阶段共 26 条数据行 × 5 个路径列（`130` 个真实文件路径；早期 M0 原始运行解释器记录保留于 `artifacts/m0/stage_report.md`）。[复算指南](reproduce.md)

3. **阶段哈希闭环校验**：
   ```bash
   sha256sum -c artifacts/m8/output_sha256.txt
   ```
   根级入口 `STATUS.md`、`README.md` 与 `SHA256SUMS.txt` 已完整同步（`artifacts/m5/` 与 `configs/project.json` 保持 M5 预解锁锁定哈希 0 字节改动）。
