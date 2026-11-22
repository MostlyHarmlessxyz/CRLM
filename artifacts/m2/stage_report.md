# M2 · GSE245552 全发现队列（V1.1 审计整改版）

2026-11-21 重建验收。39 个公开 GSM 全部导入并完成统一质量控制、官方 HGNC 规范化基因轴富集与两阶段层级细胞标注：16 位患者共有 190,190 个输入细胞，主规则保留 **173,450** 个细胞，剔除 **16,740** 个低质量/极端多联体/高线粒体细胞（剔除率 **8.80%**）。每个 GSM 具有独立的标本键和 barcode 前缀，全队列无重复完整源 barcode。`Pt05` 和 `Pt06` 各有两份 LM；`Pt08` 各有两份 PT 与 LM，均严格绑定各自 `lesion_uid` 与 `donor_uid`。[M0 临床与身份来源](../m0/manifest_review.md)

## 患者与配对病灶覆盖（宽谱系与核心细类）

| 层级 | 群体 | PT–LM 两侧 ≥20 | 两侧 ≥50（主要） | 两侧 ≥100 | NC–NL 两侧 ≥50 | 四组织 Q 共同 ≥50 |
|---|---|---:|---:|---:|---:|---:|
| `lineage` | `T/NK` | 13 | **12** | 12 | 1 | 1 |
| `lineage` | `Mono/Mac/DC` | 12 | **11** | 9 | 1 | 1 |
| `lineage` | `B/Plasma` | 8 | **8** | 5 | 1 | 0 |
| `subtype` | `CD8 T cell` | 12 | **10** | 6 | 1 | 1 |
| `subtype` | `CD4 T cell` | 12 | **12** | 9 | 1 | 1 |
| `subtype` | `monocyte/macrophage` | 12 | **11** | 9 | 1 | 1 |
| `myeloid_fine` | `Macrophage` | 11 | **11** | 7 | 0 | 0 |
| `subtype` | `plasma cell` | 8 | **6** | 5 | 1 | 0 |
| `subtype` | `B cell` | 3 | **1** | 0 | 0 | 0 |

解耦生理性共表达 `CSF3R + S100A8/A9` 但缺乏 `FCGR3B/CEACAM8` 的单核/巨噬细胞后，`Mono/Mac/DC` 主门槛（`>=50` 细胞/病灶）配对患者由 8 对恢复至 **11 对**（纯 `Macrophage` 同样达 **11 对**），且即使在对齐 `GSE164522` 的 `n_genes_visible >= 800` 敏感性门槛下仍保持 **10 对**（`T/NK = 12 对`、`B/Plasma = 8 对`）。`Pt18` 四组织在 `T/NK` 与 `Mono/Mac/DC` 上均满足 `>=50` 细胞支持；`Pt18` 的 LM 仅有 7 个 `B/Plasma` 细胞，故四组织 `B/Plasma` 有效 N 为 0。[逐患者资格及缺失原因](donor_tissue_lineage_coverage.csv) · [逐对比患者名单](paired_support.csv) · [QC 敏感性](qc_sensitivity.csv)

## 数据质量控制、两阶段层级注释与双免疫分母

| 组织 | 输入细胞 | 保留细胞 | 剔除合计 | `<300g/500c` | 谱系内双低 | 极端多联体 | `pct_mt > 20%` | 全已解析免疫分母 | 核心三谱系分母 | Nonimmune（保留） | Mixed（保留） | Unknown（保留） |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PT | 77,494 | 72,670 | 4,824 | 4,520 | 205 | 99 | 0 | 42,977 | 38,630 | 15,914 | 11,997 | 1,782 |
| LM | 83,090 | 77,571 | 5,519 | 4,982 | 477 | 60 | 0 | 35,190 | 31,998 | 27,998 | 12,490 | 1,893 |
| NC | 21,879 | 16,476 | 5,403 | 890 | 33 | 0 | 4,480 | 9,916 | 9,658 | 4,338 | 1,776 | 446 |
| NL | 7,727 | 6,733 | 994 | 689 | 298 | 7 | 0 | 4,248 | 4,089 | 475 | 1,808 | 202 |
| **合计** | **190,190** | **173,450** | **16,740** | **11,081** | **1,013** | **166** | **4,480** | **92,331** | **84,375** | **48,725** | **28,071** | **4,323** |

各库基因上限 `high_genes` 严格截断于本库可测基因数内（最高 `29,922.0`），并在宽谱系内部计算 `log1p` 离散度，彻底消除双峰混合文库（如淋巴细胞占优的 `GSM7844820`，上皮细胞保留 `464 / 477 = 97.3%`）对正常上皮细胞的误伤，同时剔除 166 个极端高计数/高基因数多联体。33 个肿瘤 GSM 缺失 `MT-` 特征（`pct_mt` 严格保留 `NaN`），`Pt18` NL（`GSM7844850`，仅 `1/7,442` 细胞检出 `MT-`）自动标记为 `visible_axis_depleted_source_anomaly`，其余 5 个具备生理性 `MT-` 分布的邻近库（4 个 NC + `GSM7844843` NL）严格剔除 4,480 个 `pct_mt > 20%` 受损细胞。[各库 QC 与排除](qc_and_exclusions.csv) · [分布图](qc_distributions.png)

在 173,450 个保留细胞中，已解析免疫细胞共 **92,331** 个（核心三谱系 **84,375** 个）：`T/NK` **46,328** 个（`CD8 T cell` 14,503、`CD4 T cell` 20,511、`Double_neg_or_unsplit_T` 10,039、`NK cell` 1,275）、`Mono/Mac/DC` **14,907** 个（`monocyte/macrophage` 14,195［含 `Macrophage` 12,291、`Monocyte` 1,904］、`conventional DC` 383、`mature DC` 288、`plasmacytoid DC` 41）、`B/Plasma` **23,140** 个（`plasma cell` 18,705、`B cell` 4,435）、`Neutrophil` 6,159 个、`Mast` 1,797 个；已解析谱系的 `unresolved` 细类率为 **0.00%**。与赋类面板完全零交集的预留独立标记在本类 vs 其他已解析类的检出率分别为：`T/NK` `91.1%` vs `8.6%`、`Mono/Mac/DC` `99.3%` vs `31.8%`、`B/Plasma` `96.1%` vs `8.0%`、`Mast` `99.1%` vs `6.5%`、`Nonimmune` `97.6%` vs `19.7%`。[注释来源与双分母规则](cell_annotation_provenance.md) · [独立标记逐库检查](annotation_marker_check.csv)

## 真实输出与验收

- [39 个 H5AD](../../data/processed/m2/) 保留各库原始整数 `X`、官方 HGNC 规范化 `var` 列、两阶段层级 `lineage_m2` / `subtype_m2` / `myeloid_fine`、主规则 `qc_retained` 及 5 组敏感性布尔掩码。4 库复用经源核对的 M1 对象，其余 35 库逐库对回 16 个源坐标。[源文件](source_inputs.csv) · [获取记录](source_fetch.csv) · [对象 SHA 与资源](object_manifest.csv)
- [三层次伪批量索引](pseudobulk_manifest.csv) 含 1,131 行（每库 29 行：8 宽谱系 + 19 层级细类 + 2 髓系细类 `Macrophage`/`Monocyte`），链接每个 GSM 的 CSR `int64` 计数和共同 29,569 个核基因顺序，并附带两套标准化免疫组成分母与敏感性计数。[基因范围与读取方法](gene_universe_audit.md)
- 独立重开全部 39 个对象，核验 190,190 个唯一细胞 ID 与 173,450 个保留细胞，逐行从筛选出的原细胞独立复算固定种子抽取的 31 个共同基因（共 35,061 个伪批量基因计数）与落盘 NPZ 完全一致；1,131 行细胞数、行计数总和及 261 行患者支持表全部复算通过。[验收记录](acceptance.json)

## 重建与资源

本机 CPU 逐库串行，Python 3.12.13。39 库构建时间合计 70.61 秒，整队列汇总另用 26.54 秒，进程峰值 RSS 为 940.13 MiB。原件约 1.07 GB，39 个 H5AD 共 1.96 GB，伪批量矩阵及基因顺序约 44.15 MB。[资源记录](resource_profile.json) · [依赖锁](python-packages.txt) · [执行入口](../../scripts/process_gse245552.py)

```bash
uv venv --python 3.12.13 .venv-m2
uv pip install --python .venv-m2/bin/python -r artifacts/m2/python-packages.txt
.venv-m2/bin/python -Bu scripts/process_gse245552.py
```

本阶段报告和全部 17 项小型表/图形产物的 SHA256 见 [产物清单](output_sha256.txt)。
