# GSE245552 发现集去偏完整候选筛选草案（V1.1）

本表严格在 `GSE245552` 开发集内基于 HGNC 规范化 `29,569` 肿瘤核基因轴与检测率分箱背景扣除 `debiased_ucell_scores` 生成（主候选上限 2–3 个）。

| 谱系 | 程序 | 配对 N | 去偏均值 D_T (95% CI) | 同向率 | log1p 同向率 | 亚型内状态分量 | 纯焦点细类 (同向率) | 筛选结论 |
|---|---|---:|---|---:|---:|---:|---|---|
| Mono/Mac/DC | `Core_Lipid_SAM_R` | 11 | +0.18092 ([+0.14409, +0.22056]) | 11/11 (100.0%) | 100.0% | +0.09145 | Macrophage: 100% (N=11) | **入选主短名单** |
| T/NK | `B07_Th17_TWEAK` | 12 | -0.03275 ([-0.05209, -0.01586]) | 10/12 (83.3%) | 91.7% | -0.03118 | CD4 T cell: 83% (N=12) | **入选主短名单** |
| Mono/Mac/DC | `B19_S4_Mac_APOE` | 11 | +0.02346 ([+0.01850, +0.03032]) | 11/11 (100.0%) | 100.0% | +0.01401 | Macrophage: 100% (N=11) | 同机制次要代表 (Core_Lipid_SAM_R) |
| Mono/Mac/DC | `B19_S4_Mac_CCL18` | 11 | +0.02119 ([+0.01499, +0.02830]) | 11/11 (100.0%) | 100.0% | +0.02012 | Macrophage: 100% (N=11) | 同机制次要代表 (Core_Lipid_SAM_R) |
| Mono/Mac/DC | `B19_S4_Mac_SPP1_clean` | 11 | +0.01625 ([+0.00796, +0.02610]) | 10/11 (90.9%) | 81.8% | +0.01243 | Macrophage: 91% (N=11) | 同机制次要代表 (Core_Lipid_SAM_R) |
| Mono/Mac/DC | `B19_S4_Mono_CCL3L1` | 11 | -0.04764 ([-0.06596, -0.02765]) | 10/11 (90.9%) | 90.9% | -0.06771 | Monocyte: 100% (N=4) | focal_subtype_not_concordant |
| Mono/Mac/DC | `Resident_Kupffer` | 11 | +0.05634 ([+0.02160, +0.09606]) | 10/11 (90.9%) | 81.8% | +0.02485 | Macrophage: 64% (N=11) | focal_subtype_not_concordant | organ_resident_control_benchmark |
| B/Plasma | `HALLMARK_UNFOLDED_PROTEIN_RESPONSE` | 8 | -0.00986 ([-0.01787, -0.00119]) | 6/8 (75.0%) | 75.0% | -0.00270 | plasma cell: 62% (N=8) | composition_dominated_below_35pct_state | focal_subtype_not_concordant |
| T/NK | `CD8_Cytotoxic_Effector` | 12 | +0.03293 ([+0.00151, +0.05976]) | 9/12 (75.0%) | 83.3% | +0.01377 | CD8 T cell: 50% (N=12) | focal_subtype_not_concordant |
| T/NK | `CD8_CXCL13_Dysfunction` | 12 | +0.01336 ([-0.00854, +0.03863]) | 8/12 (66.7%) | 58.3% | +0.01078 | CD8 T cell: 42% (N=12) | direction_support_below_75pct | bootstrap_ci_crosses_zero | log1p_scale_discordant | focal_subtype_not_concordant |
| Mono/Mac/DC | `B19_S4_Mac_F13A1` | 11 | -0.00055 ([-0.00639, +0.00631]) | 7/11 (63.6%) | 72.7% | +0.00071 | Macrophage: 73% (N=11) | direction_support_below_75pct | bootstrap_ci_crosses_zero | leave_one_out_direction | n20_direction | within_subtype_state_discordant |
| Mono/Mac/DC | `B19_S4_Mac_CXCL10` | 11 | -0.00597 ([-0.01929, +0.00798]) | 7/11 (63.6%) | 63.6% | +0.01434 | Macrophage: 55% (N=11) | direction_support_below_75pct | bootstrap_ci_crosses_zero | log1p_scale_discordant | within_subtype_state_discordant | focal_subtype_not_concordant |
| B/Plasma | `B_Antigen_Presentation_TLS` | 8 | +0.08051 ([+0.01515, +0.14937]) | 5/8 (62.5%) | 62.5% | +0.01501 | B cell: 67% (N=3) | direction_support_below_75pct | log1p_scale_discordant | composition_dominated_below_35pct_state | focal_subtype_not_concordant |
| B/Plasma | `Plasma_Antibody_Secretion` | 8 | -0.04966 ([-0.10670, -0.00163]) | 5/8 (62.5%) | 87.5% | -0.00001 | plasma cell: 38% (N=8) | direction_support_below_75pct | composition_dominated_below_35pct_state | focal_subtype_not_concordant |
| T/NK | `T_Progenitor_Memory` | 12 | -0.02258 ([-0.06595, +0.01910]) | 7/12 (58.3%) | 58.3% | -0.01544 | CD8 T cell: 58% (N=12) | direction_support_below_75pct | bootstrap_ci_crosses_zero | log1p_scale_discordant | focal_subtype_not_concordant |
| T/NK | `CD4_Treg_Suppressive` | 12 | -0.01325 ([-0.04123, +0.01181]) | 7/12 (58.3%) | 50.0% | -0.00391 | CD4 T cell: 42% (N=12) | direction_support_below_75pct | bootstrap_ci_crosses_zero | log1p_scale_discordant | composition_dominated_below_35pct_state | focal_subtype_not_concordant |
| Mono/Mac/DC | `Core_Chemokine_IFN_C` | 11 | -0.00926 ([-0.05287, +0.03492]) | 6/11 (54.5%) | 54.5% | +0.04936 | Macrophage: 55% (N=11) | direction_support_below_75pct | bootstrap_ci_crosses_zero | leave_one_out_direction | log1p_scale_discordant | within_subtype_state_discordant | focal_subtype_not_concordant |
| Mono/Mac/DC | `Core_MHCII_Presentation_A` | 11 | +0.00840 ([-0.05152, +0.06484]) | 6/11 (54.5%) | 54.5% | +0.00228 | Macrophage: 45% (N=11) | direction_support_below_75pct | bootstrap_ci_crosses_zero | leave_one_out_direction | log1p_scale_discordant | min_genes_800_align164_direction | within_subtype_state_discordant | focal_subtype_not_concordant |
| T/NK | `HALLMARK_INTERFERON_GAMMA_RESPONSE` | 12 | +0.00390 ([-0.00040, +0.00852]) | 6/12 (50.0%) | 66.7% | +0.00289 | CD4 T cell: 67% (N=12) | direction_support_below_75pct | bootstrap_ci_crosses_zero | log1p_scale_discordant | focal_subtype_not_concordant |

- **无监督 Donor-Consensus NMF 审计**：共提取 `18` 个共识模块（详见 [denovo_nmf_audit.csv](denovo_nmf_audit.csv)），与固定文献集合重合的模块均按已知程序收敛印证归档，不另造新名。
- **双尺度与组成敏感性**：完整多尺度对照见 [sensitivity_summary.csv](patient_effects/sensitivity_summary.csv)，零 `unresolved` 对称细类分解见 [patient_subtype_decomposition.csv](composition_results/patient_subtype_decomposition.csv)，纯细类患者效应见 [fine_subtypes/effects_summary.csv](patient_effects/fine_subtypes/effects_summary.csv)。
