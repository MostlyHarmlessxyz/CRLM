# CRLM 项目状态 · 2026-11-22（`CRLM_V1.1_AUDIT_REMEDIATED` · 8/8 外部检验全观测闭环）

## 三版本透明归档与外部验证属性

- **当前主交付版本**：`CRLM_V1.1_AUDIT_REMEDIATED`（冻结标识 `CRLM_V1.1_245_F1`）。其中 `GSE164522` 与 `GSE203608` 因初版曾访问表达矩阵，严格界定为**已披露前序暴露的外部稳健性重评估**；**`GSE178318`（6 对同步切除 IV 期 MSS 患者）冻结前仅审阅公开论文与 GEO 样本/条形码元数据（零前序表达矩阵或效应暴露），系在 M5 密码学锁定（`5ce11557...`）及 `5b01467` 提交后首次获授权解锁原始 UMI 矩阵并零改动运行，构成前瞻零表达暴露外部验证**（跨队列未经基因型去重，供者级非重叠状态保持 `unproven`）。
- **不可变原版基线归档**：`CRLM_V1_245_F1` 完整封存于 `artifacts/archive/CRLM_V1_245_F1/`（`archive_sha256.txt`）、`models/archive/CRLM_V1_245_F1/` 与 `embeddings/archive/CRLM_V1_245_F1/`（`ARCHIVE_SHA256.txt`）。
- **纯客观 Bug 修复对照轨**：`CRLM_V1_245_F1_techfix`（含补齐 `GSE178318` 的 8/8 实测与独立子袋 XPU 编码）完整落盘于 `artifacts/techfix_v1_f1/`。

## 正式研究阶段（`M0–M8` 全链路完成，3 个外部队列共 8/8 预登记检验零缺失）

| 阶段 | 状态 | 核心实证结果与可复核产物 |
|---|---|---|
| M0 环境与输入 | completed | 39 个 GSM、16 位患者、13 个肿瘤配对；本机隔离复核 106 项合成检查通过；含官方 HGNC 快照、`scGPT` 语料重叠与外部暴露主日志。[M0 报告](artifacts/m0/stage_report.md) |
| M1 最小真实导入 | completed | `Pt18` 四组织（22,179 细胞）稀疏导入、源坐标抽查与两阶段无冲突层级注释核验通过。[M1 报告](artifacts/m1/stage_report.md) |
| M2 全开发集 | completed | 39 库 `190,190` 输入细胞，库自适应 QC 保留 **`173,450`** 个细胞（已解析免疫谱系内部 `unresolved = 0.00%`）；`>=50` 细胞门槛下 `T/NK`、`Mono/Mac/DC`（及纯 `Macrophage`）、`B/Plasma` 配对患者分别为 **12、11、8 对**。[M2 验收](artifacts/m2/stage_report.md) |
| M3 患者程序与计数 | completed | 20 桶检测率匹配去偏 `UCell` 评估全部 19 个完整程序（`19/19` 对照零碰撞）+ 供者共识 cNMF + 有限样本 Student-$t$ / $2^N$ 符号翻转配对负二项：筛出 `Mono/Mac/DC::Core_Lipid_SAM_R`（`11/11` 升，`+0.18092, p=0.00098`）与 `T/NK::B07_Th17_TWEAK`（`10/12` 降，`−0.03275, p=0.03857`，状态占 `95.2%`）；证伪淘汰 `B/Plasma::UPR`（伪批量 `0` DEG）。[M3 报告](artifacts/m3/stage_report.md) |
| M4 全袋分类比较 | completed | A750 XPU（`torch 2.15.0+xpu`）72 细胞预检（最大绝对误差 `2.09×10⁻⁷`）后编码全袋 `70,628` 个肿瘤免疫细胞；修复细类饥饿后 `T/NK PCA_subtype_balanced` 恢复至 **`0.990`**（纯 `CD4 = 1.000`，纯 `CD8 = 0.974`），髓系 `program_elastic = 0.935`；完成 200 次置换零分布、`256×10` 子袋方差（`0.065–0.275`）、内层调参 `HA/CTA`、修复版 `MixMIL`（`All3Lineages = 1.000`）及 `15/15` 全集迁移头。[M4 报告](artifacts/m4/stage_report.md) |
| M5 原子冻结 | completed | 版本 `CRLM_V1.1_245_F1`（在获取 `GSE178318` 表达矩阵前密码学锁定，全程 0 字节改动）：锁定 2 个主候选、8 项预登记检验（4 个 Holm 家族）、冻结 `matched_controls.csv`、15 个迁移头及 72 项全链路输入 SHA256。[冻结记录](artifacts/m5/freeze_record.md) · [M5 报告](artifacts/m5/stage_report.md) |
| M6 外部复现与前瞻验证 | **completed（8/8 全观测闭环）** | **`T/NK::B07_Th17_TWEAK` 在 `GSE164522` 四组织观测交互 `I = −0.04918`（`1+/9−`, `N=10`, `Holm p=0.0430 <= 0.05`）显著通过校正**，纯 `CD4 T` 在 `164`（`1+/9−, −0.02938`）与**前瞻 `178`（`1+/5−, −0.02067`）**均高比例同向下调；**`Core_Lipid_SAM_R` 在前瞻 `178` 达 `5+/1−`（`+0.09936`，纯巨噬 `+0.10183` `[+0.02595, +0.19140]`, `5/1`）**，在 `164`（`5+/3−`, `I` `4+/1−`；共线 `Mac_SPP1_clean` `7+/1−`）与 `203608` 纯巨噬（`3+/1−`）方向一致（`0/2` 项主肿瘤 `D_T` 达显著）；**双外部队列（`164` 与前瞻 `178`）各 `15/15` 迁移头 `100%` HVG 覆盖零拒判**，`164` 上 `PCA` 全面不劣于或优于 `scGPT`，`178` 上 `T/NK` `PCA 0.917` 显著高于 `scGPT_mean 0.583`（`Δ=+0.333 [+0.083, +0.556]`；`10/12` 对 CI 跨 `0` 且 `B/Plasma` 注意力 `scGPT 0.667 > PCA 0.625`）。[M6 报告](artifacts/m6/stage_report.md) |
| M7 空间审计与验证方案 | 已完成空间审计及纯湿实验方案（湿实验执行数为 `0`） | `GSE335552`（3 患者 / 4 张 Xenium 切片）在 `M397T` 救回 `3,236` 个（四片合计 `3,584` 个）接触态 `CD8 T` 细胞，多尺度（`400->50 µm`）与半格平移（`p=0.001`）确认距离差衰减，`<15 µm EPCAM+` 肿瘤巢分层在 `M397T`（`+7.67` vs `+0.71 µm, p=0.093`）与 `M242T`（`+3.09` vs `+4.45 µm`）呈双向差异；完成四组织盲法 `SPP1+CXCL9` 双色原位杂交 + 人源髓系条件培养基（最低可行健康供者单核诱导巨噬 MDM / 增强版新鲜手术原代 TAM）→ 自体活体 `CD3+CD8+` T 细胞 Transwell 趋化（`NT_Ctrl / SPP1_KD / SPP1_Rescue_Endo / rhSPP1_Ext` 及 `CXCL9/10` 中和，对齐 `Sun Immunity 2026` 胞内 `SPP1–TRIM21–SOCS1` IFN-γ 负反馈；**已完成实际湿实验数为 0**）可操作方案。[M7 报告](artifacts/m7/stage_report.md) |
| M8 终稿与主图交付 | completed | 三版本并列终稿报告、4 幅确定性出版级主图（含前瞻 `GSE178318` 逐患者配对与真实三版本双队列迁移对比，PNG + PDF + 9 份逐元素全等源表）、局限偏差记录与 26 行全核实复算索引。[最终报告](artifacts/m8/final_report.md) · [M8 报告](artifacts/m8/stage_report.md) |

## 更新与治理规则

所有阶段状态严格依据真实落盘产物、运行日志与 `sha256sum -c` 门禁更新；`configs/project.json`、`artifacts/m5/`、`models/m4/full/` 与 `src/crlm/` 严格保持 `GSE178318` 解锁前的 M5 锁定哈希快照（确保 `sha256sum -c artifacts/m5/frozen_inputs_sha256.txt` 防篡改证明，冻结后完成的 M6 全闭环与 M8 状态以本文件及对应阶段报告为准）。`GSE225857` 因缺供者配对映射继续排除。
