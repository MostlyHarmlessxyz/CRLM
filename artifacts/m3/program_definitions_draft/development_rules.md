# GSE245552 开发集完整假设筛选、去偏评分与小样本配对统计规则（V1.1）

本规则对齐《CRLM 研究方案》§5.1–5.4、§7.4、§11 及 `IMPLEMENTATION_PLAN.md` M3 规范，100% 限定在唯一开发集 `GSE245552`（16 位患者、39 份样本、173,450 个经 M2 分层注释与 QC 保留的细胞）内执行。

## 1. 四大类候选程序版图（共 19 项完整程序 + 受控无监督共识 NMF）

| 免疫谱系（含分层细类） | 登记完整程序 | 原始文献依据与机制角色 |
|---|---|---|
| **Mono/Mac/DC**（细类：`Macrophage`、`Monocyte`、`cDC`、`pDC`、`mature DC`） | `Core_Lipid_SAM_R`、`B19_S4_Mac_SPP1_clean`、`B19_S4_Mac_CCL18`、`B19_S4_Mac_APOE`、`Core_Chemokine_IFN_C`、`B19_S4_Mac_CXCL10`、`Core_MHCII_Presentation_A`、`B19_S4_Mono_CCL3L1`、`B19_S4_Mac_F13A1`、`Resident_Kupffer` | 涵盖 H1 基质/脂质重塑 `R`、趋化/IFN `C`、抗原呈递 `A`、Bill *Science* 2023 Table S4 全部髓系状态块（客观剥离跨谱系 T/NK/B/上皮污染标志物，另保留未剥离 `B19_S4_Mac_SPP1` 作初版对照）及**肝驻留 Kupffer 对照**；Sathe 2023 `B06` S5/S6 另作通路交集与非癌疾病背景参照 |
| **T/NK**（细类：`CD8 T cell`、`CD4 T cell`、`DN/unsplit T`、`NK cell`） | `CD8_CXCL13_Dysfunction`、`CD8_Cytotoxic_Effector`、`T_Progenitor_Memory`、`CD4_Treg_Suppressive`、`B07_Th17_TWEAK`、`HALLMARK_INTERFERON_GAMMA_RESPONSE` | 涵盖 H2 肿瘤反应/耗竭、细胞毒效应、祖细胞/中心记忆、调节性 T 抑制程序、**245 原论文 `B07` `Th17–TWEAK` 参照**及 MSigDB IFN-γ 应答 |
| **B/Plasma**（细类：`B cell`、`plasma cell`） | `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`、`Plasma_Antibody_Secretion`、`B_Antigen_Presentation_TLS` | 涵盖内质网未折叠蛋白应激、恒定链/浆细胞分泌及滤泡 B 细胞抗原呈递/TLS 相关程序 |
| **受控无监督补充** | 三大谱系 `Donor-Consensus NMF`（`denovo_nmf_audit.csv`） | 肿瘤病灶均衡微池化 × 10 随机种子 NMF 共识聚类 + Kotliar-2019 标准化基因回归提取 top-50 基因；与固定文献集合收敛的模块记为已知程序收敛印证，不换名造新词 |

## 2. 基因轴规范化与深度去偏主评分公式

1. **HGNC 规范化共同核基因轴**：经 `src/crlm/hgnc_alias_map.py` 将原符号唯一映射至 HGNC Approved Symbol，锁定 `GSE245552` 四组织共同核基因轴 `245_tumor_canonical_nuclear`（`29,569` 个核基因，彻底排除 `NC/NL` 独有的 `353` 个线粒体/核糖体基因造成的计数偏移），并同步核验三队列共同核基因轴 `245_164_203608_shared_nuclear`（`18,733` 个基因）。
2. **主评分指标（`debiased_ucell_scores`）**：
   - 全测量轴 Mann–Whitney UCell（`max_rank = 1500`，对正表达取降序平均秩、对测量零赋精确全轴同值平均秩 $(K_c + 1 + G)/2$ 后截断至 `1501`）；
   - 在 `GSE245552` 肿瘤文库各谱系内计算基因检测率与平均表达分桶（`n_bins = 20`，`reference_gene_stats.csv`），从 `245_164_203608_shared_nuclear` 非成员池按同桶对称卡尺抽取 `n_ctrl = 10` 组背景对照基因（锁定于 `matched_controls.csv`），主分数为 $\tilde{S}_{c,P} = U'_c(P) - U'_c(P_{\text{ctrl}})$，从数学上消除 $K_c < 1500$ 的秩预算膨胀伪相关；
   - 同步输出 `ucell_raw`、`log1p_cp10k_debiased` 与初版 `legacy_positive`（`positive_midrank_1500`）双尺度及历史对照。

## 3. 五重预设定短名单筛选门禁（主候选总数上限 2–3 个）

1. **样本量门禁**：主门槛（每文库该谱系 `>= 50` 个 QC 细胞）下完整 PT–LM 配对患者数 $N \ge 8$；
2. **一致性与区间门禁**：非零患者与均差同向比例 $\ge 75\%$，患者级 10,000 次 bootstrap 95% CI 不跨 0，且留一患者（LOO）均差符号完全不变；
3. **双尺度与 QC 门槛稳健门禁**：在连续表达尺度 `log1p_cp10k_debiased` 下均差同向且非零同向率 $\ge 70\%$，同时 `n20`、`n100` 与 `min_genes_800_align164`（对齐 `GSE164522` 基因数下限）均差保持同向；
4. **非纯组成伪影门禁**：在零 `unresolved` 分层细类共同支持对称分解中，`within_subtype_state` 与总 `D_T` 同向且占比 $\ge 35\%$，且在对应纯焦点细类（如 `Macrophage`、`CD8 T cell`、`CD4 T cell`、`plasma cell`，门槛 `>= 20` 细胞）内 $N \ge 5$、均差同向且非零同向率 $\ge 70\%$（客观拦截纯属亚型比例变化的伪候选）；`Resident_Kupffer` 单列作肝驻留器官背景参照；
5. **同机制共线去重门禁**：同谱系内同机制簇或患者 $D_T$ $|\rho| \ge 0.75$ / 基因重叠系数 $\ge 0.30$ 的通过候选，按预设洁净度层级仅保留最优代表进入主检验家族（其余归入次要敏感性代表，避免同机制重复计入稀释外部 Holm 功效），跨谱系总入选上限为 3 项。

## 4. 小样本稳健配对伪批量基因推断

对三大谱系（`T/NK`、`Mono/Mac/DC`、`B/Plasma`）及四大核心细类（`CD8 T cell`、`CD4 T cell`、`Macrophage`、`plasma cell`）分别拟合 `~ donor_uid + tissue` 负二项 GLM（`poscounts` 大小因子，严格屏蔽 `_genewise_converged / _MAP_converged / _LFC_converged` 未收敛基因并移除配对设计下恒为空操作的 `cooks_filter`），以残差自由度 $df = N_{\text{pairs}} - 1$ 的双侧 Student-$t$ 分布计算主 `pvalue_t` / `padj_t`，并辅以无分布假设的 $2^{N_{\text{pairs}}}$ 全枚举双侧精确配对符号翻转置换检验（核基因 `log2(CP10k + 0.5)` 差值，`pvalue_signflip` / `padj_signflip`）交叉印证。
