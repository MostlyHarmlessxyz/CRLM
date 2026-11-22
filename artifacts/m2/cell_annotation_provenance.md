# GSE245552 发现集统一 QC、两阶段层级标注与双免疫分母（V1.1 审计整改版）

## 计数与统一 QC 纳入规则

每库沿用已核对的 GSM、作者标本号、原始 barcode 与交付整数计数。首先基于 `src/crlm/annotation_m1.py` 对所有细胞执行两阶段确定性标注，再在同一 GSM 内按宽谱系（`n >= 50` 时取该谱系内分布，否则回退至整库分布）计算可见计数和检出基因数的 `log1p` 中位数与 `1.4826 * MAD` 界限，并将基因上限截断于各库交付维度内（`high_genes <= native_genes`）。主要排除条件（`qc_retained == False`）按优先级互斥记录于 `obs.qc_exclusion`：

1. `empty`：`n_counts_visible == 0`。
2. `below_min_300g_500c`：`n_genes_visible < 300` 或 `n_counts_visible < 500` 硬下限剔除低复杂度碎片。
3. `within_lineage_low`：同时低于所属宽谱系（或整库回退）`median - 3.0 * 1.4826 * MAD` 的双指标低端细胞。
4. `joint_high_multiplet`：硬上限多联体（`n_genes_visible > 8,500` 或 `n_counts_visible > 80,000`），或同时高于谱系内 `median + 3.5 * 1.4826 * MAD` 且满足生物学多联体下界（`n_genes_visible > 6,000` 或 `n_counts_visible > 40,000`），避免在淋巴细胞占优的肿瘤库（如 `GSM7844820`）中误杀正常高表达肿瘤上皮细胞。
5. `high_pct_mt_gt20`：仅对具备真实生理性线粒体分布的邻近正常组织库（`pct_mt_scope == "visible_axis"`，共 5 个 GSM）剔除 `pct_mt > 20.0%` 的受损细胞。

肿瘤交付轴（33 个 GSM）缺失全部 13 个 `MT-` 基因，`pct_mt` 严格保留为 `NaN`（`unavailable_no_MT_features`）；`Pt18` 的邻近肝 `GSM7844850` 虽含 13 个 `MT-` 行，但 `7,441 / 7,442`（`99.987%`）细胞的 `MT-` 计数为 0，自动标记为 `visible_axis_depleted_source_anomaly`（`pct_mt_usable_for_qc = False`），不施加虚假的 `pct_mt` 截断。所有细胞及其 5 组敏感性布尔掩码（`qc_pass_g200`、`qc_pass_g300_primary`、`qc_pass_g500`、`qc_pass_g800` 对齐 `GSE164522` 的 800 基因门槛、`qc_pass_no_mt_filter`）均完整保存在 H5AD 的 `obs` 中，配对患者人数敏感性汇总于 [QC 规则敏感性](qc_sensitivity.csv)。队列 CD45 磁选及重混协议见 [M0 原始来源核查](../m0/manifest_review.md) 与 [B07 正文](../../references/papers/B07_2024.pdf)。

## 两阶段层级身份标记与显式双联体剥离

1. **第一阶段（粗谱系解耦与跨谱系双联体剥离）**：将 `CSF3R` 移出中性粒细胞锚点，保留给生理性共表达 `CSF3R + S100A8/A9` 但缺乏 `FCGR3B/CEACAM8` 的经典单核/巨噬细胞，真性中性粒细胞由 `FCGR3B/CEACAM8/CXCR2/CMTM2` 锚定；同时将显式跨谱系混合表达（`CD3+B`、`CD3+多髓系`、`B+多髓系`、`PTPRC+EPCAM/HNF4A+` 等）统一剥离至 `Mixed`（`subtype_m2 = "Cross_lineage_doublet"`），不进入任何单谱系状态或免疫组成分母。
2. **第二阶段（谱系内零 `unresolved` 层级细分）**：
   - `T/NK`：`CD8 T cell`、`CD4 T cell`、`Double_neg_or_unsplit_T`、`NK cell`（优先保护肝脏 `CD3D/E/G/TRAC == 0` 但携带无菌 `TRBC1/2` 转录本的 `NK cell`）；
   - `Mono/Mac/DC`：`monocyte/macrophage`（并在 `obs.myeloid_fine` 及 `level == "myeloid_fine"` 进一步无损细分为 `Macrophage` 与 `Monocyte`）、`conventional DC`、`plasmacytoid DC`、`mature DC`；
   - `B/Plasma`：`B cell`、`plasma cell`；
   - `Nonimmune`：`epithelial`、`hepatocyte`、`endothelial`、`fibroblast`、`mural`。

程序评分和 PT–LM 效应均未参与上述规则。另采用与赋类面板及辅助标记**完全零交集**的预留独立标记（`HELD_OUT_CORROBORATION_GENES`）逐库核验本类与其他已解析类的检出率，结果记录于 [独立标记检查表](annotation_marker_check.csv)。[标记规则](../../src/crlm/annotation_m1.py) · [M2 入口](../../scripts/process_gse245552.py)

## 估计单位与双免疫组成分母

每个 GSM 各自计数；同患者同组织的多个病灶（如 `Pt08` 的 2 PT + 2 LM、`Pt05/Pt06` 的 2 LM）保留各自 `lesion_uid`。由于肿瘤与正常文库经 CD45+ 与 CD45− 磁珠分选后重混，全捕获细胞分母受实验富集比干扰；M2 在 [伪批量索引](pseudobulk_manifest.csv) 与 [各库 QC 表](qc_and_exclusions.csv) 中同步记录两套标准化免疫组成分母：

- `n_total_resolved_immune`（`prop_total_resolved_immune`）：`T/NK + Mono/Mac/DC + B/Plasma + Neutrophil + Mast`；
- `n_core_three_lineage_immune`（`prop_core_three_lineage_immune`）：`T/NK + Mono/Mac/DC + B/Plasma`（跨队列与仅含核心三谱系的外部子集直接对齐）。

`Nonimmune`、`Mixed` 与 `Unknown` 的免疫组成比例严格记为 `NaN`。细胞内连续状态主门槛要求每个标本×谱系（或细类）至少 50 个保留细胞，并同步列出 20、100 的敏感性。[患者覆盖](donor_tissue_lineage_coverage.csv) · [配对支持](paired_support.csv)
