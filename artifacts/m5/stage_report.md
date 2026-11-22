# M5 · V1.1 方法学重构版候选程序与分类模型原子冻结报告（`CRLM_V1.1_245_F1`）

> [!IMPORTANT]
> **三版本透明并列与历史暴露声明**
> 初版基线 `CRLM_V1_245_F1`（`artifacts/archive/CRLM_V1_245_F1/`）与纯技术勘误 `CRLM_V1_245_F1_techfix`（`artifacts/techfix_v1_f1/`）已完整封存并记录于 [external_access_log.csv](external_access_log.csv)。本冻结包 `CRLM_V1.1_245_F1` 的候选筛选、20 桶 × 10 重复检测率匹配对照基因集与 15 个全开发集分类头均 100% 在唯一开发集 `GSE245552` 内重做，外部评估严格界定为**已披露历史暴露的独立队列复现与方法学稳健性验证**。

## 1. 冻结主候选摘要（来自 `GSE245552` 五重门禁筛选）

| 主候选 (`lineage::program`) | 冻结预期方向 | 开发集同向患者 | `mean_D_T` [95% CI] | 双侧符号 `p` | 焦点纯细类 `mean_D_T` (同向率) | 亚型内状态占比 |
|---|---|---:|---:|---:|---:|---:|
| `Mono/Mac/DC::Core_Lipid_SAM_R` | 上升 (LM > PT, +1) | 11/11 (100.0%) | +0.18092 [+0.14409, +0.22056] | 0.00098 | `Macrophage`: +0.16857 (100.0%) | 50.5% |
| `T/NK::B07_Th17_TWEAK` | 下降 (LM < PT, -1) | 10/12 (83.3%) | -0.03275 [-0.05209, -0.01586] | 0.03857 | `CD4 T cell`: -0.03888 (83.3%) | 95.2% |

- **同机制共线次要候选与器官驻留基准（不计入主检验多重比较家族）**：`Mono/Mac/DC::B19_S4_Mac_APOE` (secondary_collinear_with:Core_Lipid_SAM_R)、`Mono/Mac/DC::B19_S4_Mac_CCL18` (secondary_collinear_with:Core_Lipid_SAM_R)、`Mono/Mac/DC::B19_S4_Mac_SPP1_clean` (secondary_collinear_with:Core_Lipid_SAM_R)、`Mono/Mac/DC::Resident_Kupffer` (organ_resident_control)（详见 [secondary_sensitivity_candidates.csv](secondary_sensitivity_candidates.csv)）。
- **初版 `B/Plasma::HALLMARK_UNFOLDED_PROTEIN_RESPONSE` 淘汰依据**：去偏后开发集 `6/8 = 75.0%`（`p = 0.289`），`72.6%` 总差由 `13/13` 患者浆细胞比例暴跌驱动，纯 `plasma cell` 内仅 `5/8 = 62.5%`（`-0.00308`），未通过 Gate 4 组成分解与焦点细类门禁。

## 2. 预登记检验矩阵与全开发集迁移分类模型

- **预登记检验（共 `8` 项，分属 4 个独立 Holm 校正家族，每族 `m = 2`，`alpha = 0.05`）**：
  - `PRIMARY_DT_164_V1_1`（`GSE164522` 配对肿瘤差 `D_T = LM - PT`）
  - `BACKGROUND_I_164_V1_1`（`GSE164522` 四组织合格子集 `Q` 交互 `I = D_T_Q - D_N_Q`）
  - `SECONDARY_DT_178_V1_1`（`GSE178318` 次要方向证据；使用条款待核实前保留 `raw_p = NA`、Holm 向量置 `1.0`）
  - `SECONDARY_DT_203608_V1_1`（`GSE203608` 次要方向与精度证据）
- **冻结全开发集迁移分类头（`15/15` 项全部验收通过）**：`T/NK`、`Mono/Mac/DC`、`B/Plasma` × `{program_elastic, PCA_mean, PCA_attention, scGPT_mean, scGPT_attention}`，全部基于 `245_164_shared_nuclear`（`18,815` 个共享核编码基因，在 `GSE164522`、`GSE203608` 与 `scGPT` 词表上实现 `512/512 = 100.0%` PCA HVG 覆盖，彻底消除初版 4 个模型拒判）。

## 3. 冻结清单与防篡改校验

- **时间与版本**：`2026-11-21T17:49:04+00:00` · `CRLM_V1.1_245_F1` · Git HEAD `91a589f5ec58736753e5fd71138e0733e3bc8ab1`
- **配置自洽性**：`configs/project.json` 与 [development_config.json](development_config.json) 在冻结时刻原子同步（`research_version = "V1.1_AUDIT_REMEDIATED"`，`freeze_version = "CRLM_V1.1_245_F1"`）。
- **核心文件**：[冻结记录](freeze_record.md) · [主要候选](frozen_candidate_registry.csv) · [预登记检验](registered_tests.csv) · [验证规则](validation_rules.md) · [匹配对照集](matched_controls.csv) · [冻结模型](frozen_models.csv) · [外部暴露日志](external_access_log.csv) · [全链路输入哈希](frozen_inputs_sha256.txt) · [阶段产物哈希](output_sha256.txt)
