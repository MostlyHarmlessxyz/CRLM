# GSE245552 测量范围与原计数

导入器逐库将实际 feature 符号与 M0 审计的组织基因集合比对，检查源 feature ID、符号和 barcode 的唯一性。保存时 `X` 为细胞 × 原始交付基因的 CSR int32 计数，`var` 同时保留原 ID 和符号。细胞 ID 为 GSM 加作者 barcode。每库输出先重读核对稀疏类型、矩阵轴、样本映射和总计数；Pt18 复用已经对回源矩阵的 M1 对象，其余样本每库抽查 16 个源坐标。SHA、GEO URL 和各库对象见 [源文件清单](source_inputs.csv) 及 [输出清单](object_manifest.csv)。

PT、LM 各有 29,569 个基因，NC、NL 各有 29,922 个，四组织共同可见轴为 29,569。邻近组织额外的 353 个符号含全部 13 个 `MT-` 基因。缺失行保留为不可测，肿瘤 `pct_mt` 为缺失。各库 `n_counts_visible` 与 `n_genes_visible` 均针对本库交付轴。逐基因归属可查 [共同集合](../../starter/audit_data/GSE245552_feature_intersection.txt) 和 [邻近组织特有轴](../../starter/audit_data/GSE245552_normal_only_features_353.txt)。

## 伪批量

每个 GSM 的宽谱系与其二级子组单独将通过 QC 的细胞原计数相加，随后保留上述共同可测基因并按 [共同基因顺序](../../data/processed/m2/pseudobulk/gene_symbols.txt) 排列。[伪批量索引](pseudobulk_manifest.csv) 的 `count_matrix` 与零基 `matrix_row` 可直接索引 SciPy CSR NPZ，`n_cells` 和 `common_axis_visible_count_sum` 分别是细胞数和这一行计数和。宽谱系行与二级行各自构成完整分区，统计合并时只选一个层次。一侧某亚型为零细胞时，索引保留该行并标为不满足状态覆盖。

该分支支持原计数的患者感知分析。对于多个病灶，合并原计数表示捕获细胞混合，病灶等权的连续程序均值使用细胞对象另算。跨事件、治疗与独立取材时序依据及限制继承自 [M0 临床来源](../m0/manifest_review.md)。
