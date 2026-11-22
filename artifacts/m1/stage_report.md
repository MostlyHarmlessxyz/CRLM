# M1 · Pt18 四组织稀疏导入（V1.1 审计整改版）

2026-11-21 重建验收。GSE245552 的 Pt18 同时满足 PT–LM 配对和四组织要求，输入为 1 位患者、4 个标本、4 个公开 GSM、22,179 个细胞。每个 GSM 导出一个可重读的 H5AD。四个对象共有 40,958,372 个非零整数计数条目。本轮升级集成了官方 HGNC 规范化基因轴映射、两阶段层级细胞注释与显式双联体剥离，以及邻近组织线粒体源异常自动识别。

| 组织 | GSM | 输入与保留细胞 | 基因数 | 可见计数中位数/细胞 | 检出基因中位数/细胞 | 检出 MT 计数的细胞 | `pct_mt_scope` | 可用于 M2 MT 过滤 |
|---|---|---:|---:|---:|---:|---:|---|---|
| 原发灶 PT | GSM7844847 | 4,639 | 29,569 | 2,675 | 1,309 | 不可测 | `unavailable_no_MT_features` | 否 |
| 肝转移灶 LM | GSM7844848 | 7,009 | 29,569 | 2,297 | 1,442 | 不可测 | `unavailable_no_MT_features` | 否 |
| 邻近结直肠 NC | GSM7844849 | 3,089 | 29,922 | 4,322 | 1,352 | 885 (28.65%, 684 >20%) | `visible_axis` | 是 |
| 邻近肝 NL | GSM7844850 | 7,442 | 29,922 | 3,189 | 1,734.5 | 1 (0.013%) | `visible_axis_depleted_source_anomaly` | 否 |

## 输入、HGNC 规范轴与 QC 能力审计

[选样表](sample_selection.csv) 取自 [M0 样本表](../../data/manifests/m0/sample_manifest_validated.csv)。每个 GSM 的矩阵、feature 和 barcode 都来自 GEO，12 件原文件的 SHA 与逐 GSM 来源记录一致。[输入路径、URL 与 SHA](source_inputs.csv)

[四个 H5AD](../../data/processed/m1/) 在 `X` 保存细胞 × 基因的 CSR `int32` 计数，`var` 保留源 `source_feature_id` 与原始 `gene_symbol`，并新增基于本地 SHA 锁定的官方 HGNC 完整表映射列（`canonical_symbol`、`hgnc_id`、`hgnc_match_type`、`in_excluded_353`、`in_tumor_canonical_nuclear_245`、`in_shared_nuclear_245_164`、`in_shared_nuclear_245_164_203608`），不改变原始矩阵维度或顺序。重读验证了形状、轴顺序、身份、层级亚型、HGNC 规范符号、稀疏格式、非零数、可见计数总和及每库 16 个源坐标（共 64 个）。[逐库矩阵验收与 H5AD SHA](matrix_audit.csv)

输入为作者交付、经过基因预筛的整数计数。PT、LM 共享 29,569 个核基因轴，NC、NL 共享 29,922 个基因轴；邻近组织多出的 353 个基因包含全部 13 个 `MT-` 基因，故 PT、LM 全部细胞的 `pct_mt` 严格保留为 `NaN`（绝不伪造为 0）。在邻近组织中，NC（`GSM7844849`）呈现生理性线粒体分布（885 个细胞检出 `MT-` 计数，684 个细胞 `pct_mt > 20%`，最高 `49.96%`，标记为 `visible_axis`），而 NL（`GSM7844850`）虽在 `features.tsv.gz` 含有 13 个 `MT-` 行，但 `7,441 / 7,442`（`99.987%`）细胞的 `MT-` 计数为 0（仅 1 个细胞检出 148 UMI，`pct_mt = 22.19%`），由导入器自动识别为 `visible_axis_depleted_source_anomaly`（`pct_mt_usable_for_qc = False`），避免下游把源端近乎清零的 `MT-` 轴误当作已完成正常线粒体过滤。

10 组文献种子同时报告原符号精确覆盖与 HGNC 规范化后覆盖：例如 `B19_S4_Mac_CXCL10` 在 PT/LM 由原符号 `90/94` 提升至规范化 `93/94`（NC/NL 达 `94/94`，找回 `C19orf66→SHFL`、`C9orf66→DOCK8-AS1`、`RARRES3→PLAAT4`；仅 `RPS6KC1` 属 353 个仅正常组织可测基因），`B19_S4_Mac_SPP1` 在 PT/LM 由 `160/164` 提升至 `162/164`（找回 `FAM129B→NIBAN2`、`FAM213B→PRXL2B`；`LONP1` 属 353 个仅正常组织可测基因，`ST20` 不在交付轴）。[逐基因范围](gene_coverage.csv) · [种子基因覆盖](seed_feature_coverage.csv)

M1 仅剔除总计数为零的细胞（本批为 0），保留分布第 1、第 99 百分位尾部复查标记与完整捕获细胞分母；全队列硬下限、谱系内离散度异常值、极端多联体及可用 `visible_axis` 线粒体过滤在 M2 执行。[QC 分布、尾部数与 MT 可测性](qc_capabilities.csv) · [协议和临床来源](../m0/manifest_review.md)

## 两阶段层级谱系覆盖与双免疫分母

| 组织 | T/NK | Mono/Mac/DC | B/Plasma | Neutrophil | Mast | Nonimmune | Mixed | Unknown | 全已解析免疫分母 | 核心三谱系分母 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PT | 434 | 525 | 156 | 149 | 15 | 2,242 | 540 | 578 | 1,279 | 1,115 |
| LM | 262 | 112 | 7 | 3 | 3 | 5,897 | 525 | 200 | 387 | 381 |
| NC | 668 | 62 | 842 | 1 | 49 | 1,013 | 397 | 57 | 1,622 | 1,572 |
| NL | 3,204 | 763 | 317 | 290 | 5 | 490 | 1,998 | 375 | 4,579 | 4,284 |

两阶段层级注释器（`src/crlm/annotation_m1.py`）先解耦生理性共表达 `CSF3R + S100A8/A9` 但缺乏 `FCGR3B/CEACAM8` 的经典单核/巨噬细胞与真性中性粒细胞，并将显式跨谱系双联体（`CD3+B`、`CD3+多髓系`、`B+多髓系`、`PTPRC+EPCAM/HNF4A+` 等）剥离至 `Mixed`（`Cross_lineage_doublet`）；随后在各已解析谱系内做确定性层级细分（`T/NK`：`CD8 T cell`、`CD4 T cell`、`Double_neg_or_unsplit_T`、`NK cell`；`Mono/Mac/DC`：`monocyte/macrophage`［并细分 `Macrophage`/`Monocyte`］、`conventional DC`、`plasmacytoid DC`、`mature DC`；`B/Plasma`：`B cell`、`plasma cell`），已解析谱系的 `unresolved` 比例为 `0.00%`。覆盖表同时记录全捕获保留分母、全已解析免疫分母（含 `Neutrophil` 与 `Mast`）及核心三谱系免疫分母（`T/NK + Mono/Mac/DC + B/Plasma`）。每库每大类抽取至多 3 个细胞，94 个抽查细胞的 925 条原始标记计数均与 H5AD 逐项一致。[层级注释代码](../../src/crlm/annotation_m1.py) · [逐库覆盖与双免疫分母](lineage_coverage.csv) · [抽查细胞](marker_spotchecks.csv)

例如，PT 的 `GSM7844847:R5218-3_AGACGTTTCTTCATGT-1` 检出 `C1QA=1`、`CD14=1`、`FCN1=1`、`VCAN=8`、`CD68=1`、`CSF3R=1` 且 `FCGR3B=0`、`CEACAM8=0`，由旧规则下的假性 `Mixed` 恢复为 `Mono/Mac/DC`（`monocyte/macrophage` / `Monocyte`）；NL 的 `Mono/Mac/DC` 覆盖由 586 升至 763（含 607 个巨噬细胞、53 个单核细胞、85 个 cDC、9 个 pDC、9 个成熟 DC），肝 `NK cell` 恢复至 337 个。

## 执行与资源

本机 i5-13400，CPU 逐库串行，Python 3.12.13。12 件压缩原文件共 119.03 MiB。原件已就位时，导入、HGNC 映射、两阶段注释和重读校验共 10.232 秒，进程峰值 RSS 816.96 MiB。4 个 H5AD 共 216.62 MiB，压缩 MTX 流式读取，未生成解压矩阵临时文件。完整环境与代码/词表哈希见 [依赖锁](python-packages.txt) 和 [资源记录](resource_profile.json)，逐库日志见 [import.log](import.log)。

从项目根目录重建（随机抽查种子 `20261121`）：

```bash
uv venv --python 3.12.13 .venv-m1
uv pip install --python .venv-m1/bin/python -r artifacts/m1/python-packages.txt
.venv-m1/bin/python -Bu scripts/import_gse245552.py \
  --manifest data/manifests/m0/sample_manifest_validated.csv --donor GSE245552:Pt18
```

[导入脚本](../../scripts/import_gse245552.py) 会校验并复用本地一致的 GEO 原件与官方 HGNC 表。[M1 报表与日志 SHA](output_sha256.txt)

## M2 交接

其余 35 个 GSM 从同一份 M0 样本表处理，复用 `import_gse245552.py` 的稀疏读取、HGNC 规范轴富集、`pct_mt_scope` 异常识别与两阶段层级注释，并在 M2 施加全队列统一 QC（`n_genes_visible >= 300`、`n_counts_visible >= 500`、谱系内 `log1p` MAD 下界、极端多联体上限剔除，以及 `visible_axis` 正常组织 `pct_mt <= 20%` 过滤）与敏感性梯度核查。
