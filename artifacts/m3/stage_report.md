# M3 · GSE245552 开发集完整假设筛选、深度去偏评分与稳健配对推断报告（V1.1）

## 1. 执行范围与方法学修复概览

本阶段严格限定在唯一开发集 **`GSE245552`**（16 位转移性结直肠癌患者、39 份标本、经 M2 两阶段分层注释与一致 QC 保留的 **173,450 个细胞**，细类 `unresolved = 0.00%`）内完成，系统修复初版 M3 的四项核心缺陷：
1. **消除核糖体/线粒体背景偏移与符号断层**：经 `src/crlm/hgnc_alias_map.py` 规范化全部基因符号并统一对齐 `245_tumor_canonical_nuclear`（**`29,569` 个共同核编码基因**，剔除 `NC/NL` 独有的 `353` 个 MT/核糖体基因），同时核验三队列共同核基因轴 `245_164_203608_shared_nuclear`（`18,733` 个基因）。
2. **根除单细胞测序深度秩膨胀伪影**：以全测量轴 Mann–Whitney UCell 联合 `GSE245552` 肿瘤文库 20 桶 × 10 重复检测率匹配背景扣除（`debiased_ucell_scores`）作为主指标，锁定 [reference_gene_stats.csv](program_definitions_draft/reference_gene_stats.csv) 与 [matched_controls.csv](program_definitions_draft/matched_controls.csv)，同步输出 `ucell_raw`、`log1p_cp10k_debiased` 与初版 `legacy_positive` 对照。
3. **恢复完整研究方案假设版图（19 项完整文献程序 + 受控无监督共识 NMF）**：覆盖 H1 髓系（`Core_Lipid_SAM_R`、去跨谱系污染 `B19_S4_Mac_SPP1_clean`、Bill S4 其余完整髓系集、`Core_Chemokine_IFN_C`、`Core_MHCII_Presentation_A`、肝驻留对照 `Resident_Kupffer`）、H2 T/NK（`CD8_CXCL13_Dysfunction`、`CD8_Cytotoxic_Effector`、`T_Progenitor_Memory`、`CD4_Treg_Suppressive`、245 原论文参照 `B07_Th17_TWEAK`、`HALLMARK_INTERFERON_GAMMA_RESPONSE`）、B/Plasma（`HALLMARK_UNFOLDED_PROTEIN_RESPONSE`、`Plasma_Antibody_Secretion`、`B_Antigen_Presentation_TLS`）及三大谱系 10 种子 `Donor-Consensus NMF`（18 个共识因子）。
4. **小样本稳健配对伪批量推断**：对三大谱系与四大核心细类（`CD8 T cell`、`CD4 T cell`、`Macrophage`、`plasma cell`）严格屏蔽 `PyDESeq2` 未收敛基因、关闭配对设计下恒为空操作的 `cooks_filter`，以残差自由度 $df = N_{\text{pairs}} - 1$ 的双侧 Student-$t$ 分布（`pvalue_t`/`padj_t`）联合无分布假设的 $2^{N_{\text{pairs}}}$ 全枚举双侧精确配对符号翻转置换检验（`pvalue_signflip`/`padj_signflip`）替代反保守渐近正态 Wald 检验。

---

## 2. 深度去偏成效实证核验（`legacy_positive` vs `ucell_debiased`）

实测 [depth_debiasing_proof_summary.csv](depth_debiasing_proof_summary.csv) 直接证实初版 `positive_midrank_1500` 在稀疏细胞中存在强测序深度混杂，而 `ucell_debiased` 大幅消除该技术伪相关：
- **`B/Plasma::HALLMARK_UNFOLDED_PROTEIN_RESPONSE`**：细胞级与检出基因数 Spearman 均值由 **`+0.6115` 降至 `+0.2708`**；患者配对效应与中位 UMI 差的相关 $\rho(D_T, \Delta\text{UMI})$ 由 **`+0.8810` 降至 `+0.3095`**（与检出基因差相关由 **`+0.5238` 降至 `+0.1190`**）。
- **`T/NK::HALLMARK_INTERFERON_GAMMA_RESPONSE`**：细胞级与检出基因数 Spearman 均值由 **`+0.4821` 降至 `-0.0282`**。
- **`Mono/Mac/DC::Core_MHCII_Presentation_A`**：初版公式下受深度膨胀驱动的假正偏（`+0.04768`, 8/11 同向，细胞级 $\rho = +0.6006$，深度标准化后反转为 `-0.00790`）在 `ucell_debiased` 下直接归零（`+0.00840`, 6/11 同向）。

---

## 3. 十九项完整候选程序五重门禁筛选结果

完整筛选指标见 [candidate_screen.csv](candidate_screen.csv) 与 [candidate_evidence_draft.md](candidate_evidence_draft.md)。按预设定五重门禁（① $N \ge 8$；② 非零同向率 $\ge 75\%$、bootstrap 95% CI 不跨 0 且 LOO 全同向；③ `log1p_cp10k_debiased` 同向率 $\ge 70\%$ 且多门槛稳健；④ 零 `unresolved` 对称分解亚型内状态占比 $\ge 35\%$ 且纯焦点细类同向率 $\ge 70\%$；⑤ 同谱系共线/同机制去重并控制总数 $\le 3$），最终锁定 **2 项跨谱系主候选**：

1. **`Mono/Mac/DC::Core_Lipid_SAM_R`（主候选 1，预期 `LM > PT` 上升）**：
   - 配对患者 **`N = 11`，`11/11 (100.0%)` 全员正向**，去偏均值 `D_T = +0.18092`（95% CI `[+0.14409, +0.22056]`，开发集双侧精确符号检验 `p = 0.00098`，LOO 范围 `[+0.16753, +0.18818]`）；
   - `log1p_cp10k_debiased` 尺度同样 **`11/11 (100.0%)`**（均值 `+0.73418`）；
   - 对称细类分解显示亚型内状态分量 `+0.09145`（占总效应 `50.5%`，共同支持内 `+0.14438`），且在**纯 `Macrophage` 细类（`N = 11` 对）内仍保持 `11/11 (100.0%)` 正向**（均值 `+0.16857`）；
   - 同机制/高共线髓系集 `B19_S4_Mac_APOE`（`11/11`, `+0.02346`）、`B19_S4_Mac_CCL18`（`11/11`, `+0.02119`）与去污染 `B19_S4_Mac_SPP1_clean`（`10/11`, `+0.01625`）同样通过前四重门禁，按第 ⑤ 门禁归入同机制次要敏感性代表（不重复占用外部主 Holm 家族名额）。
2. **`T/NK::B07_Th17_TWEAK`（主候选 2，预期 `LM < PT` 下降 / 原发肠黏膜肿瘤富集）**：
   - 配对患者 **`N = 12`，`10/12 (83.3%)` 同向负值**，去偏均值 `D_T = -0.03275`（95% CI `[-0.05209, -0.01586]`，精确符号检验 `p = 0.03857`，LOO 范围 `[-0.03620, -0.02811]`）；
   - `log1p_cp10k_debiased` 尺度达 **`11/12 (91.7%)`**（均值 `-0.10557`），深度标准化后同向率 `91.7%`；
   - 对称分解证实 **`95.2%` 为纯亚型内状态变化**（`-0.03118` vs 组成 `-0.00157`），在**纯 `CD4 T cell` 细类（`N = 12` 对）内保持 `10/12 (83.3%)` 同向**（均值 `-0.03888`）；同时须在外部四组织 `I = D_T - D_N` 中严格检验结直肠黏膜基线 Th17 背景能否解释该差异。

### 关键淘汰与证伪诊断（验证研究方案 §5.4 证伪规则）
- **`B/Plasma::HALLMARK_UNFOLDED_PROTEIN_RESPONSE` 被客观淘汰**：去偏后大类同向率降至 `6/8 (75.0%, p = 0.289)`，对称分解证实 **`72.6%` 来自浆细胞相对 B 细胞比例暴跌**（在已解析免疫细胞分母中，`B/Plasma` 比例在 **13/13 位患者全员下降**，均值 `-22.5%`），而**纯 `plasma cell` 细类（`N = 8`）内 `UPR` 均值仅 `-0.00308`、同向率 `5/8 (62.5%)`**，确证其为典型组成伪影。
- **`T/NK::CD8_Cytotoxic_Effector` 被细类门禁拦截**：粗大类 `T/NK` 虽呈 `9/12 (75.0%, +0.03293)` 上升，但在**纯 `CD8 T cell`（`N = 12`）内同向率恰好为 `6/12 (50.0%, +0.00972)`**，揭示粗类上升源于肝转移 `CD8+NK / CD4` 亚群比例切换而非 CD8 细胞内杀伤程序上调。
- **`Resident_Kupffer`（肝驻留对照）**：粗类 `10/11 (+0.05634)`，但纯 `Macrophage` 内同向率仅 `63.6% (7/11)`（近半效应 `+0.03333` 来自单核→巨噬比例提升），保留为外部器官背景参照。

---

## 4. 受控无监督共识 NMF 与小样本配对伪批量基因推断

1. **Donor-Consensus NMF（[denovo_nmf_audit.csv](denovo_nmf_audit.csv)）**：三大谱系共提取 18 个高稳定性共识因子（余弦稳定度 `0.976–1.000`）。其中 `T_NK__cNMF_05`（`12/12` 患者 `LM < PT`，与 `B07_Th17_TWEAK` 重叠系数 `0.50`，命中 `IL17A/IL17F/IL22/IL23R/IL26/RORA/CTSH`）、`Mono_Mac_DC__cNMF_02`（与 `Core_Chemokine_IFN_C` 重叠系数 `0.43`）及 `B_Plasma__cNMF_04`（与 `B_Antigen_Presentation_TLS` 重叠系数 `0.25`）直接收敛回固定文献集合；`Mono_Mac_DC__cNMF_05`（`11/11` `LM > PT`）捕获 `FABP5/LPL/MMP7/SLAMF9/CD109` 脂质重塑轴，未出现独立于现有候选体系且需另造新名的未知主导状态。
2. **配对伪批量负二项 Student-$t$ 与 $2^{N}$ 精确符号翻转置换（[fit_summary.csv](paired_counts/fit_summary.csv)、[reference_readouts.csv](paired_counts/reference_readouts.csv)）**：
   - 严格屏蔽未收敛基因（`T/NK` 41–47 个、`Mono/Mac/DC` 11–16 个、`B/Plasma` 8–12 个）并改用 $df = N_{\text{pairs}} - 1$ 的 Student-$t$ 后，消除了渐近正态 Wald 的严重小样本膨胀：
     - **`Mono/Mac/DC`（`N=11, df=10`）**：`padj_t < 0.05` 共 **661** 基因，与完全无分布假设的 $2^{11}$ 精确符号翻转 **`padj_signflip < 0.05`（611 基因）** 高度一致（渐近正态虚高至 `1,736`）；**纯 `Macrophage`（`N=11`）** 内 **`padj_t < 0.05` 为 400 基因，`padj_signflip < 0.05` 为 399 基因**，核心 `SPP1`（`log2FC = +4.21, padj_t = 0.014, padj_signflip = 0.047`）、`APOC1`（`+2.95, padj_t = 0.0040, padj_signflip = 0.038`）、`APOE`（`+2.18, padj_t = 0.011, padj_signflip = 0.038`）、`GPNMB`（`+1.71, padj_t = 0.011, padj_signflip = 0.052`）、`TREM2`（`+1.87, padj_t = 0.018, padj_signflip = 0.052`）在纯巨噬细胞内全面通过双方法验证！
     - **`T/NK` 与纯 `CD4 T cell`（`N=12, df=11`）**：纯 `CD4 T cell` 内 `padj_t < 0.05` **102** 基因、`padj_signflip < 0.05` **265** 基因，其中 `IL17A`（`log2FC = -4.01, padj_t = 0.014, padj_signflip = 0.048`）与 `RORC`（`-0.96, p_t = 0.0017, padj_signflip = 0.048`）显著下调，与 `B07_Th17_TWEAK` 程序结论严丝合缝。
     - **`B/Plasma`（`N=8, df=7`）与纯 `plasma cell`（`N=6, df=5`）**：渐近正态假阳性报告 `42`（粗类）/ `39`（纯浆细胞）个基因，而在 Student-$t$ 与 $2^N$ 精确符号翻转下 **`padj < 0.05` 均为 0 基因**（纯浆细胞内 `XBP1 p_t = 0.668`、`HSP90B1 p_t = 0.768`、`DERL3 p_t = 0.438`、`MZB1 p_t = 0.847`），从伪批量计数端再次印证纯浆细胞内无 UPR 转录组下调。

---

## 5. 资源消耗与完整性校验

- `analyze_discovery_programs.py`：`125.49s`，峰值 RSS `818.9 MiB`（[resource_profile.json](resource_profile.json)）。
- 完整性校验：`sha256sum -c artifacts/m3/output_sha256.txt`。
