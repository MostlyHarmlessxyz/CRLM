# M6 · 外部队列程序复现、四组织器官背景与三版本对照（V1.1 `CRLM_V1.1_245_F1`）

验证版本为 **`CRLM_V1.1_245_F1`**（基于披露 `CRLM_V1_245_F1` 历史暴露后的开发集方法学去偏重构）。8 项预登记外部检验分布在 4 个独立 Holm 家族（每族 `m = 2`），其中 6 项在 `GSE164522` 与 `GSE203608` 获得实测估计，2 项（`GSE178318`）因作者数据使用条款待书面确认按冻结规则保留 `raw_p = NA`（仅 Holm 校正输入置 `1.0`）。[登记检验结果](registered_tests/registered_test_results.csv) · [家族汇总](registered_tests/registered_families.csv) · [三版本并列对照](three_version_evidence_comparison.md)

## 一、GSE164522 主外部队列与四组织背景（156,511 洁净细胞，18,815 共享核基因）

| 检验 ID | 对比 | 谱系 :: 程序 | 患者 N | 同向 / 反向 | 平均差 [95% bootstrap CI] | 双侧 sign p | 族内 Holm p |
|---|---|---|---:|---:|---|---:|---:|
| `PRIMARY_DT_164_V1_1_C1` | `D_T` | `Mono/Mac/DC::Core_Lipid_SAM_R` | 8 | 5 升 / 3 降 | `+0.01677` `[−0.02096, +0.05813]` | 0.7266 | 0.7266 |
| `PRIMARY_DT_164_V1_1_C2` | `D_T` | `T/NK::B07_Th17_TWEAK` | 10 | 7 降 / 3 升 | `−0.00923` `[−0.02178, +0.00366]` | 0.3438 | 0.6875 |
| `BACKGROUND_I_164_V1_1_C1` | `I` | `Mono/Mac/DC::Core_Lipid_SAM_R` | 5 | 4 正 / 1 负 | `+0.07425` `[+0.00009, +0.15752]` | 0.3750 | 0.3750 |
| `BACKGROUND_I_164_V1_1_C2` | `I` | `T/NK::B07_Th17_TWEAK` | 10 | **9 负 / 1 正** | **`−0.04918`** `[−0.07201, −0.02687]` | **0.0215** | **0.0430** |

1. **四组织背景交互检验通过（`BACKGROUND_I_164_V1_1_C2`）**：`T/NK::B07_Th17_TWEAK` 在正常肝高于正常结肠（`D_N_Q = +0.03995`，`10/10` 全正），而在肝转移灶相对原发灶逆转下调（`D_T_Q = −0.00923`；等价地 `PT − NC = +0.03404`、`LM − NL = −0.01514`），使同一患者四组织观测交互差 `I = D_T_Q − D_N_Q = −0.04918`（9/10 为负，`Holm p = 0.0430 <= 0.05`）在全部双尺度（`ucell_debiased` 与 `log1p_cp10k_debiased`）及统一规则重注释（`rule_annotation_m1`）分支下稳定显著。
2. **谱系异质性与焦点细类验证**：宽 `T/NK` 的 `D_T` 受肝内 `CD8 T cell` 正向偏移（`+0.02003`, `10/0`）部分对冲，而在生物学对应的纯 **`CD4 T cell`** 内，`B07_Th17_TWEAK` 的 `D_T` 达 **`−0.02938`（9 降 / 1 升，`raw p = 0.0215`）**；在洁净作者标注纯 `Macrophage` 内 `Core_Lipid_SAM_R` 的 `D_T` 为 `+0.02847`（`5/3`；`rule_annotation_m1` 分支为 `−0.00589`, `3/5`），共线次要集 `B19_S4_Mac_SPP1_clean` 达 **`7 升 / 1 降`（`+0.00936`, `raw p = 0.0703`；四组织 `I` 为 `5/0`, `p = 0.0625`）**。
3. **V1 历史伪影证伪**：作为未进入主族的探索性证伪对照，`B/Plasma::UPR` 在洁净四组织患者（`N=9`）的 `D_T_Q`（`−0.01381`）与 `D_N_Q`（`−0.01291`）近似抵消（`I = −0.00090`, `95% CI [−0.00721, +0.00473]`, `6 正 / 3 负`），结合对称分解（`V1.1` 下组成分量占 `77.3%`）证实其跨部位差异主要源于器官背景与浆细胞比例变化。[164 详报](GSE164522/stage_report.md) · [焦点细类总表](focal_subtype_effects_all_cohorts.csv) · [统一注释对照](uniform_annotation_effects_summary.csv)

## 二、次要外部来源（GSE203608 与 GSE178318）及迁移分类衔接

- **`GSE203608`（47,525 洁净体内细胞，18,733 共享核基因）**：宽 `Mono/Mac/DC::Core_Lipid_SAM_R` 均值为 `+0.01713`（`2/2`, `N=4`, `direction_matches = False`, `p_H = 1.0`；剔除 Stage III `Patient 49` 后为 `+0.07141`, `2/1`），在纯 **`Macrophage`** 内（洁净作者标注与 `rule_annotation_m1` 双分支一致）均为 **`3 升 / 1 降`（`+0.03051` / `+0.01134`）**；`T/NK::B07_Th17_TWEAK` 在 `>=50` 门槛下仅 `N=2`（`−0.00850`, `1/1`），在 `rule_annotation_m1` `>=20` 门槛下为 `0 升 / 3 降`（`−0.02169`）。[UPenn 详报](GSE203608/stage_report.md)
- **`GSE178318`（`terms_pending`）**：未读取任何表达值，两项次要检验保留 `raw_p = NA`、`holm_adjusted_p = 1.0`。[条款状态](GSE178318/access_status.md)
- **外部迁移分类器（15/15 全部激活，6 个 PCA 头 512/512 = 100.0% HVG 覆盖，0 特征拒判）**：基于同一套洁净 `GSE164522` 肿瘤标本（21 个 PT/LM 标本、57 个 `>=50` 细胞的合格全袋、`68,050` 个经 Intel Arc A750 XPU 编码的肿瘤免疫细胞，`285/315` 个袋×方法产出预测、其余 `30` 个因 `<50` 细胞未入袋）与 `18,815` 共同核基因轴执行全袋推理：
  - `T/NK`（10 供者 / 20 部位 / 10 对）：`PCA_attention = 0.980` `[0.920, 1.000]`（10/10 对）、`PCA_mean = 0.970` `[0.900, 1.000]`（10/10 对），点估计均高于 `scGPT_attention = 0.880` 与 `scGPT_mean = 0.810`（配对池化差 `PCA_mean - scGPT_mean = +0.160 [+0.020, +0.250]`，`PCA_attention - scGPT_attention = +0.100 [0.000, +0.270]`），`program_elastic = 0.500`；
  - `Mono/Mac/DC`（9 供者 / 17 部位 / 8 对）：`PCA_attention = 0.889` `[0.630, 1.000]`、`scGPT_attention = 0.806`、`PCA_mean = 0.750`、`scGPT_mean = 0.625`（`PCA_attention - scGPT_mean = +0.264 [+0.049, +0.444]`）、`program_elastic = 0.583`；
  - `B/Plasma`（10 供者 / 19 部位 / 9 对）：`program_elastic` 彻底消除初版符号反转并升至 **`0.856`** `[0.675, 1.000]`（8/9 对 `LM > PT`），`PCA_attention = 0.844` `[0.683, 0.971]`、`PCA_mean = 0.667` 均高于 `scGPT_mean = 0.578`（`PCA_attention - scGPT_mean = +0.267 [+0.040, +0.517]`）与 `scGPT_attention = 0.333`。[模型推理详报](GSE164522/model_predictions/stage_report.md) · [跨来源证据卡](cross_cohort_evidence_cards.md) · [访问审计](external_access_log.csv)
