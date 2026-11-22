# GSE245552 测量范围、HGNC 规范轴与分层伪批量（V1.1 审计整改版）

导入器逐库将实际 feature 符号与 M0 审计的组织基因集合及本地 SHA256 锁定的官方 HGNC 完整表比对，检查源 `source_feature_id`、原始 `gene_symbol`、规范化 `canonical_symbol` 和 barcode 的唯一性。保存时 `X` 为细胞 × 原始交付基因的 CSR `int32` 计数；`var` 保留原始 `source_feature_id` 与 `gene_symbol`，并附带确定性 HGNC 规范化列（`canonical_symbol`、`hgnc_id`、`hgnc_match_type`、`in_excluded_353`、`in_tumor_canonical_nuclear_245`、`in_shared_nuclear_245_164`、`in_shared_nuclear_245_164_203608`）。`Pt18` 四组织复用已对回源矩阵的 M1 对象，其余 35 个样本每库抽查 16 个源坐标。SHA、GEO URL 和各库对象见 [源文件清单](source_inputs.csv) 及 [输出清单](object_manifest.csv)。

PT、LM 各有 29,569 个核基因，NC、NL 各有 29,922 个基因，四组织共同可见核轴为 29,569（HGNC 规范化后与 `GSE164522` 共享核基因为 `18,815`，三队列 `245∩164∩203608` 共享核基因为 `18,733`）。邻近组织额外的 353 个符号（含全部 13 个 `MT-` 基因及若干核糖体/特定非编码基因）保留在 `excluded_353_canonical` 清单中，跨组织程序评分与迁移模型严格排除这 353 个基因；肿瘤 `pct_mt` 保持 `NaN`，绝不填 0。[共同集合](../../starter/audit_data/GSE245552_feature_intersection.txt) · [邻近组织特有轴](../../starter/audit_data/GSE245552_normal_only_features_353.txt) · [HGNC 规范轴清单](../../data/manifests/hgnc/PROVENANCE.json)

## 三层次伪批量矩阵

每个 GSM 按以下三个层次将 `qc_retained == True` 的细胞原计数相加，保留共同可测的 29,569 个核基因并按 [共同基因顺序](../../data/processed/m2/pseudobulk/gene_symbols.txt) 排列（每库 29 行，39 库共 1,131 行）：

1. `level == "lineage"`（第 `0..7` 行，共 8 行）：宽谱系完整互斥分区（`T/NK`、`Mono/Mac/DC`、`B/Plasma`、`Neutrophil`、`Mast`、`Nonimmune`、`Mixed`、`Unknown`）；
2. `level == "subtype"`（第 `8..26` 行，共 19 行）：两阶段层级细类完整互斥分区（含 `CD8 T cell`、`CD4 T cell`、`Double_neg_or_unsplit_T`、`NK cell`、`monocyte/macrophage`、`conventional DC`、`plasmacytoid DC`、`mature DC`、`B cell`、`plasma cell` 等，零 `unresolved`，且每个父谱系下细类计数之和严格等于父谱系行）；
3. `level == "myeloid_fine"`（第 `27..28` 行，共 2 行）：`monocyte/macrophage` 进一步无损细分的 `Macrophage` 与 `Monocyte`（两者之和严格等于 `subtype == "monocyte/macrophage"` 行）。

[伪批量索引](pseudobulk_manifest.csv) 的 `count_matrix` 与零基 `matrix_row` 可直接索引 SciPy CSR NPZ，供 M3 大谱系与核心细类（`CD8 T cell`、`CD4 T cell`、`Macrophage`、`plasma cell`、`B cell`）配对原计数建模及有效患者数核查直接调用。
