# CRLM 项目状态 · 2026-11-22（`CRLM_V1.1_AUDIT_REMEDIATED`）

## 三版本透明归档与当前版本属性

- **当前主交付版本**：`CRLM_V1.1_AUDIT_REMEDIATED`（冻结标识 `CRLM_V1.1_245_F1`，完整修复初版 7 项方法学缺陷后的去偏重评估版本；鉴于初版曾访问外部队列，其外部检验明确定位为已披露历史暴露的独立队列稳健性重评估）。
- **不可变原版基线归档**：`CRLM_V1_245_F1` 的 `M3–M8` 阶段产物与根级快照、256 细胞抽样嵌入及模型已完整封存于 `artifacts/archive/CRLM_V1_245_F1/`（`archive_sha256.txt`）、`models/archive/CRLM_V1_245_F1/` 与 `embeddings/archive/CRLM_V1_245_F1/`（`ARCHIVE_SHA256.txt`）。
- **纯客观 Bug 修复对照轨**：`CRLM_V1_245_F1_techfix`（锁定初版候选、方向、`positive_midrank_1500` 公式及初版模型权重，仅修复 HGNC 符号、共享核空间重标定、双细胞/死细胞污染及 `MixMIL` `/sqrt(Q)` 与初始化截断 Bug）完整落盘于 `artifacts/techfix_v1_f1/`。

## 正式研究阶段（`M0–M8` 全链路完成）

| 阶段 | 状态 | 核心实证结果与可复核产物 |
|---|---|---|
| M0 环境与输入 | completed | 39 个 GSM、16 位患者、13 个肿瘤配对；本机隔离复核 106 项合成检查通过；补充 `scGPT` 语料重叠与官方 HGNC 快照审计。[M0 报告](artifacts/m0/stage_report.md) |
| M1 最小真实导入 | completed | `Pt18` 四组织（22,179 细胞）稀疏矩阵导入、源坐标抽查与两阶段无冲突层级注释核验通过。[M1 报告](artifacts/m1/stage_report.md) |
| M2 全开发集 | completed | 39 库 `190,190` 输入细胞，库自适应 QC 保留 **`173,450`** 个细胞（已解析谱系 `unresolved = 0.00%`）；`>=50` 细胞门槛下 `T/NK`、`Mono/Mac/DC`（及纯 `Macrophage`）、`B/Plasma` 配对患者分别为 **12、11、8 对**；生成三层次 `int64` 伪批量。[M2 验收](artifacts/m2/stage_report.md) |
| M3 患者程序与计数 | completed | 20 桶检测率匹配去偏 `UCell` 评估全部 19 个完整程序 + 供者共识 cNMF + 有限样本 Student-$t$ 配对负二项与 `log₂(CP10k+0.5)` 精确符号翻转：筛出 `Mono/Mac/DC::Core_Lipid_SAM_R`（`11/11` 升，`+0.18092, p=0.00098`）与 `T/NK::B07_Th17_TWEAK`（`10/12` 降，`−0.03275, p=0.03857`，状态占 `95.2%`）；证伪淘汰组成驱动的 `B/Plasma::UPR`（伪批量 `0` DEG）。[M3 报告](artifacts/m3/stage_report.md) |
| M4 全袋分类比较 | completed | A750 XPU 72 细胞预检（最大绝对误差 `2.09×10⁻⁷`）后编码全袋 `70,628` 个肿瘤免疫细胞；5 折供者分层 CV 下 `T/NK PCA_subtype_balanced`（`CD8/CD4` 等权质心，`16p/28s/12` 对）达 **`0.990`**（纯 `CD4 = 1.000`，纯 `CD8 = 0.974`），髓系 `program_elastic = 0.935`；完成 200 次置换零分布、`256×10` 子袋方差、内层调参 `HA/CTA`、修复版 `MixMIL`（`All3Lineages = 1.000`）及 `15/15` 全集迁移头。[M4 报告](artifacts/m4/stage_report.md) |
| M5 原子冻结 | completed | 版本 `CRLM_V1.1_245_F1`：锁定 2 个主候选、8 项预登记检验（4 个 Holm 家族）、冻结 `matched_controls.csv`、15 个迁移头及 72 项全链路输入 SHA256，与 `configs/project.json` 保持同哈希快照并透明披露前序暴露。[冻结记录](artifacts/m5/freeze_record.md) · [M5 报告](artifacts/m5/stage_report.md) |
| M6 外部复现与迁移 | 已完成可用来源，`178` 待条款 | **`T/NK::B07_Th17_TWEAK` 在 `GSE164522` 四组织观测交互 `I = −0.04918`（`1+/9−`, `N=10`, `raw p=0.0215`, `Holm p=0.0430 <= 0.05`）显著通过家族校正**，焦点纯 `CD4 T cell` `D_T = −0.02938`（`1+/9−`）；`Core_Lipid_SAM_R` 在 `164`（`5+/3−`, `I` `4+/1−`；共线 `Mac_SPP1_clean` `7+/1−` `D_T`, `5+/0−` `I`）与 `203608` 纯巨噬（`3+/1−`，宽髓系 `2+/2−`）未达主族显著；**`15/15` 外部迁移头全部激活（6 个 PCA 头 `512/512 = 100%` HVG 覆盖，0 特征拒判）**，`PCA` 点估计在三谱系均高于对应 `scGPT`（`T/NK 0.980/0.970`, 髓系 `0.889`, `B/Plasma program_elastic 0.856`）。`GSE178318` 保持零表达暴露 `terms_pending`。[M6 报告](artifacts/m6/stage_report.md) |
| M7 空间与功能设计 | 已完成空间审计及验证预案 | `GSE335552`（3 患者 / 4 切片）四片合计救回 `3,584` 个接触态 `CD8 T` 细胞（其中 `M397T` 救回 `3,236` 个），多尺度（`400->50 µm`）与 `<15 µm EPCAM+` 距离分层揭示 `M397T` `200 µm` 距离差集中于近上皮带（`+7.67 µm` vs 间质 `+0.71 µm, p=0.093`）；建立四组织盲法 `SPP1+CXCL9` 双染分区与同供者髓系条件液–`CD8+ T` 定向趋化（区分内源恢复 vs 生理外源 `rhSPP1`，参考 `Sun Immunity 2026`）操作预案。[M7 报告](artifacts/m7/stage_report.md) |
| M8 终稿与主图交付 | completed | 三版本并列终稿报告、4 幅主图（含 `GSE178318` 未暴露占位面板，PNG + PDF + 7 份源表）、局限与偏差记录及 24 行路径全核实的复算索引，阶段 SHA256 全绿。[最终报告](artifacts/m8/final_report.md) · [M8 报告](artifacts/m8/stage_report.md) |

## 更新与治理规则

所有阶段状态严格依据真实落盘产物、运行日志与 `sha256sum -c` 门禁更新；`configs/project.json` 与 `artifacts/m5/development_config.json` 保持 M5 预验证冻结时刻的同哈希锁定快照（确保 `sha256sum -c artifacts/m5/frozen_inputs_sha256.txt` 全链路防篡改闭环，冻结后完成的 M6/M8 最终验收状态以本文件及对应阶段报告为准）。完整保留未显著、不可测（`GSE178318` 待条款、`GSE225857` 排除）及前序外部暴露记录。
