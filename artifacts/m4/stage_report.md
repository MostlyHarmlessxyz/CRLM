# M4 · 开发集全袋细胞表示编码、去混杂消融与组织分类基准报告（Track 2 · `CRLM_V1.1`）

> [!IMPORTANT]
> 本报告对应 `CRLM_V1.1_AUDIT_REMEDIATED` 在唯一开发集 `GSE245552` 上的重评估结果（初版 `CRLM_V1_245_F1` 的 256 细胞抽样产物与模型已完整封存于 `artifacts/archive/CRLM_V1_245_F1/artifacts/m4/`、`models/archive/CRLM_V1_245_F1/` 与 `embeddings/archive/CRLM_V1_245_F1/`）。分类器作为评估转录状态可分性与器官背景印记的支撑工具，跨队列生物学复现结论由 M5–M6 冻结患者配对检验决定。

## 1. 核心实证发现与初审疑点闭环解答

1. **`T/NK` 细类平衡骤降（初版 `1.000 -> 0.521`）机制彻底查明**：
   - 初版骤降完全源于粗糙两分（`T cell` vs 稀有 `NK`）叠加 `min(n_T, n_NK) >= 3` 随机丢弃 90% 以上细胞造成的严重细胞饥饿（中位仅剩 6 个细胞）；
   - 在 `V1.1` 零 `unresolved` 细类下，保留全部 `CD8 T cell` 与 `CD4 T cell` 求质心并严格 `50/50` 等权的 **`PCA_subtype_balanced` 在全部 16 位患者 / 28 个部位上 AUROC 达 `0.990`**（12/12 配对全员正向）；**纯 `CD8 T cell`（`PCA_pure_primary`）AUROC 为 `0.974`**，**纯 `CD4 T cell`（`PCA_pure_secondary`）AUROC 为 `1.000`**；而仅用细类比例的 `composition_elastic` 为 `0.508`、仅用测序深度的 `QC_elastic` 为 `0.405`；
   - 这确证 `T/NK` 的高区分度并非亚型比例偏移、测序深度混杂或测试泄漏，而是真实的**纯 T 细胞亚型内器官微环境转录差异**（含黏膜 Th17-TWEAK vs 肝脏驻留/代谢印记；折内剔除 `NL/NC` 非免疫实质高表达基因后 `PCA_remove_organ_markers` 仍达 `1.000`，表明该信号内生于淋巴细胞转录组而非游离肝细胞 RNA 污染）。
2. **`program_elastic` 符号反转修复与深度去偏固定程序基线**：
   - 引入类别均衡患者部位权重（`sum(w_PT) == sum(w_LM) == 0.5`，消除初版未配对 PT 患者双倍权重失真）+ 配对排序损失 + 折内锁定 20 桶检测率匹配对照的 `debiased_ucell_scores` 后，髓系 `program_elastic` AUROC 由 `0.828` 提升至 **`0.935`**（11/11 配对全员正向），`T/NK` 升至 `0.595`，`B/Plasma` 彻底消除初版 `0.175` 反向失真并回归随机水平 `0.567`（与 M3 纯浆细胞无显著转移重编程结论完全一致）。
3. **全袋保留 vs `256 x 10` 子袋抽样方差审计（`subbag_256x10_sensitivity.csv`）**：
   - 全量保留 `33` 个合格肿瘤文库的 `88` 个单谱系袋（`70,388` 个细胞）与 `33` 个跨三谱系联合袋（`70,628` 个细胞）；
   - 10 次独立 `256` 细胞抽样实验证明：均值池化对子袋截断稳健（SD `0.002–0.030`），但**注意力池化受 256 细胞单次截断扰动极大**——`PCA_attention` 在 10 次抽样间的 AUROC 极差在 `B/Plasma` 达 **`0.217`**（`0.600–0.817`）、`T/NK` 达 **`0.154`**（`0.826–0.979`）、髓系达 **`0.118`**（`0.864–0.982`），`scGPT_attention` 在 `B/Plasma` 极差达 **`0.233`**，定量证实全袋评估的必要性。
4. **200 次患者级随机置换零分布（`PCA/permutation_null_summary.csv`）**：
   - 修复初版单次哈希仅翻转 3/13 对患者的偏倚后，`200` 次独立患者级 Rademacher 符号翻转零分布均值严格居中：`T/NK` 为 **`0.523 ± 0.125`**（95% 区间 `[0.287, 0.774]`，经验单侧 $p = 0.00498$），`Mono/Mac/DC` 为 **`0.515 ± 0.125`**（`[0.278, 0.746]`，$p = 0.00498$），`B/Plasma` 为 **`0.506 ± 0.140`**（`[0.225, 0.767]`，$p = 0.0199$）。
5. **复杂层级/变分 MIL 与简约均值基线对比**：
   - 在单谱系与 `All3Lineages` 全标本任务上，修复后的 `PCA32_MixMIL`（`All3Lineages = 1.000`、髓系 `0.994`、`T/NK = 0.985`）与 `PCA_mean`（`All3Lineages = 0.986`、髓系 `0.988`、`T/NK = 0.995`）并列最优，内层调参 `HA`/`CTA` 无跨谱系增量（详见 [直接前作报告](direct_prior/stage_report.md)）。

## 2. 主种子（`20261121`）折外患者组织级 AUROC 汇总表

| 模型 / 对照分支 | `T/NK` (16人/28部位/12对) | `Mono/Mac/DC` (15人/26部位/11对) | `B/Plasma` (15人/23部位/8对) | `All3Lineages` (16人/29部位/13对) |
|---|---:|---:|---:|---:|
| `QC_elastic`（可见 UMI/基因数分位数对照） | 0.405 | 0.757 | 0.633 | — |
| `composition_elastic`（全标本细类比例基线） | 0.508 | 0.870 | 0.617 | — |
| `program_elastic`（折内去偏固定程序基线） | 0.595 | 0.935 | 0.567 | — |
| `PCA_mean`（全袋 32D 均值池化） | **0.995** | 0.988 | 0.808 | 0.986 |
| `PCA_attention`（透明线性 attention） | 0.903 | 0.982 | 0.658 | 0.986 |
| `PCA_subtype_balanced`（主导双细类 50/50 等权质心） | **0.990** | 0.968 | 0.705 | — |
| `PCA_pure_primary`（纯 `CD8 T` / `Macrophage` / `plasma`） | 0.974 | **1.000** | **0.942** | — |
| `PCA_pure_secondary`（纯 `CD4 T` / `Monocyte` / `B cell`） | **1.000** | 0.968 | 0.509 | — |
| `PCA_remove_organ_markers`（折内剔除 `NL/NC` 实质基因） | **1.000** | 0.988 | 0.783 | — |
| `PCA_permuted`（200 次患者置换均值 `[2.5%, 97.5%]`） | 0.523 `[0.29, 0.77]` | 0.515 `[0.28, 0.75]` | 0.506 `[0.22, 0.77]` | — |
| `PCA32_MixMIL`（修复版全协方差变分 GLMM） | 0.985 | **0.994** | 0.783 | **1.000** |
| `scGPT_mean`（原生 512D 全袋均值池化） | 0.892 | 0.893 | **0.817** | 0.981 |
| `scGPT_attention`（原生 512D 线性 attention） | 0.644 | 0.911 | 0.300 | 0.910 |
| `scGPT32_mean`（折内 32D 投影均值池化） | 0.908 | 0.953 | **0.833** | 0.971 |
| `scGPT32_MixMIL`（折内 32D 修复版变分 GLMM） | 0.887 | 0.941 | 0.725 | 0.957 |

- [完整患者级指标](patient_metrics.csv) · [模型配对 AUROC 差值区间](paired_model_comparisons/auc_differences.csv) · [256×10 子袋敏感性](subbag_256x10_sensitivity.csv) · [折内实质掩码审计](PCA/organ_marker_mask_audit.csv) · [200 次置换零分布](PCA/permutation_null_summary.csv)

## 3. Intel Arc A750 XPU 全袋编码、归因稳定性与全集迁移头验收

1. **72 细胞 CPU vs A750 XPU 精度预检与全量缓存（`resource_profiles/`）**：
   - 在 `GSM7844847` 三谱系共 **`72` 个真实细胞**上，A750 XPU 相对 CPU 的最小余弦相似度为 **`0.99999982`**（均值 `1.0000000`），**最大逐维绝对误差仅 `2.09 × 10⁻⁷`**（均值 `1.15 × 10⁻⁸`）；跨批大小重排最小余弦为 **`0.99999988`**、最大绝对误差 **`8.94 × 10⁻⁸`**，峰值显存分配 `1.86 GiB`（[scgpt_profile.json](resource_profiles/scgpt_profile.json)）；
   - `245_164_shared_nuclear`（`18,815` 个核基因）经 HGNC 单射桥接命中 `scGPT vocab.json` 中 **`18,717` 个无碰撞词表 token**，覆盖全部非零计数的 **`99.923%`**；
   - 全部 `33` 个肿瘤文库 `70,628` 个合格免疫细胞在 A750 XPU 上耗时 `925.96s`（约 `76.3 cells/s`），缓存总大小 `135.75 MB`（`38,898` 个细胞超 `1,199` 表达基因上限并采用确定性逐细胞种子抽样），100% 通过维度、有限值与单位范数校验（[scgpt_cache_manifest.csv](resource_profiles/scgpt_cache_manifest.csv)）；预训练语料重叠核查见 `artifacts/m0/scgpt_pretraining_overlap_audit.md`（`GSE245552` / `GSE164522` / `GSE203608` 均晚于或不在 `scGPT` whole-human 2023 CELLxGENE 语料中）。
2. **折外跨种子细胞归因稳定性（`attribution_stability.csv`）**：
   - 在 `338` 组折外跨种子（`20261121` vs `20261122/20261123`）全袋比较中，`PCA_attention` 带符号贡献 Spearman 秩相关中位数为 **`1.000`**（前 10% 绝对贡献细胞集 Jaccard 中位数 **`1.000`**），`scGPT_attention` 中位数提升至 **`0.768`**（Jaccard 中位数 `0.398`，较初版 256 细胞抽样下的 `0.594 / 0.250` 明显改善但仍具高维多模态敏感性）。
3. **全开发集 15 个外部迁移头（`full_data_models.csv`）**：
   - 共登记 `490` 次折内/前作运行（`327` 次 L-BFGS-B 收敛、`120` 次内层 CV 调参完成、`40` 次 MixMIL 预定步数完成、`3` 次高维 512D 无正则投影 `scGPT` 辅助种子触达 500 步上限并诚实保留状态），共生成 `2,881` 条折外预测；
   - 在全开发集重拟合的三谱系 **15/15 个预登记迁移分类头全部收敛通过**，所有 `pca.npz` 严格锁定在 `245_164_shared_nuclear` 上（确保外部 `GSE164522` 评估 `512/512 = 100.0%` 特征覆盖、零拒判），且每个谱系目录均同步锁定冻结 `program_defs.json` 与 `matched_controls.csv`；
   - `artifacts/m4/` 根级与全部子目录文件均已写入 [output_sha256.txt](output_sha256.txt) 并通过 `sha256sum -c --quiet` 全量核验。
