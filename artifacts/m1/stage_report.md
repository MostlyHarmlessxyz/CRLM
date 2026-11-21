# M1 · Pt18 四组织稀疏导入

2026-11-21 完成。GSE245552 的 Pt18 同时满足 PT–LM 配对和四组织要求，输入为 1 位患者、4 个标本、4 个公开 GSM、22,179 个细胞。每个 GSM 导出一个可重读的 H5AD。四个对象共有 40,958,372 个非零整数计数条目。

| 组织 | GSM | 输入与保留细胞 | 基因数 | 可见计数中位数/细胞 | 检出基因中位数/细胞 | 检出 MT 计数的细胞 |
|---|---|---:|---:|---:|---:|---:|
| 原发灶 PT | GSM7844847 | 4,639 | 29,569 | 2,675 | 1,309 | 不可测 |
| 肝转移灶 LM | GSM7844848 | 7,009 | 29,569 | 2,297 | 1,442 | 不可测 |
| 邻近结直肠 NC | GSM7844849 | 3,089 | 29,922 | 4,322 | 1,352 | 885 |
| 邻近肝 NL | GSM7844850 | 7,442 | 29,922 | 3,189 | 1,734.5 | 1 |

## 输入、重读与 QC

[选样表](sample_selection.csv) 取自 [M0 样本表](../../data/manifests/m0/sample_manifest_validated.csv)。每个 GSM 的矩阵、feature 和 barcode 都来自 GEO，12 件原文件的 SHA 与逐 GSM 来源记录一致。[输入路径、URL 与 SHA](source_inputs.csv)

[四个 H5AD](../../data/processed/m1/) 在 `X` 保存细胞 × 基因的 CSR int32 计数，`var` 保留源 feature ID 与原符号，本批两者逐行相同。四库的单元键均为 GSM 加原始 barcode，`obs` 保留患者、标本、GSM、捕获标识与组织。重读验证了形状、轴顺序、身份、稀疏格式、非零数、可见计数总和及每库 16 个源坐标，共 64 个。[逐库矩阵验收与 H5AD SHA](matrix_audit.csv)

输入为作者交付、经过基因筛选的整数计数。两种肿瘤组织的测量轴相同，邻近组织多出的 353 个基因含全部 13 个 MT 基因，故 PT、LM 全部细胞的 `pct_mt` 为缺失。NC、NL 按各自可见计数为分母，MT 比例中位数均为 0，实际检出细胞数见上表。基因名按原符号精确匹配。10 组文献种子逐组织列出可测与缺失符号，仅记录基因是否在交付轴上。[逐基因范围](gene_coverage.csv) · [种子基因覆盖](seed_feature_coverage.csv)

M1 仅剔除总计数为零的细胞，本批为 0。各库可见总计数与检出基因数的第 1、第 99 百分位用于标出分布尾部以供复查，所有非空细胞计入覆盖分母。标本捕获经过作者报告的 CD45 阳性与阴性磁选及 2∶1 重混，逐库实际回收比例未知。[QC 分布、尾部数与 MT 可测性](qc_capabilities.csv) · [协议和临床来源](../m0/manifest_review.md)

## 初级谱系覆盖

| 组织 | T/NK | Mono/Mac/DC | B/Plasma | Neutrophil | Mast | Nonimmune | Mixed | Unknown |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| PT | 469 | 409 | 172 | 204 | 15 | 2,355 | 462 | 553 |
| LM | 321 | 95 | 10 | 4 | 5 | 5,955 | 414 | 205 |
| NC | 695 | 56 | 898 | 2 | 52 | 1,127 | 195 | 64 |
| NL | 3,615 | 586 | 380 | 312 | 5 | 574 | 1,526 | 444 |

标签由 15 组身份标记产生。命中一组需检出至少两个不同基因，其中至少一个为该组身份标记。同一大类内的命中合并，跨大类命中记为 Mixed，无合格组合记为 Unknown。S3 原表的基因与作者分组用于交叉参照。每库每类抽取至多 3 个细胞，保留原始计数和匹配组合。95 个细胞的 823 条标记计数均与 H5AD 一致。[初级注释代码](../../src/crlm/annotation_m1.py) · [逐库覆盖](lineage_coverage.csv) · [抽查细胞](marker_spotchecks.csv)

例如，PT 条形码 `R5218-3_CCAGCGAAGTCTCCTC-1` 检出 CD3E=6、CD3D=3、TRAC=2，支持 T/NK 初判。PT 的 `R5218-3_CGACTTCCATCCAACA-1` 同时满足浆细胞与上皮组合，保留为 Mixed。NL 有 1,526 个混合标记细胞，LM 只得到 10 个 B/Plasma 初判。M2 应结合其他文库与更细致的注释核对这些类别，再确定可用于患者配对分析的细胞支持。此处的计数只描述该患者捕获样本中的初级覆盖。

## 执行与资源

本机 i5-13400，CPU 逐库串行，Python 3.12.13。12 件压缩原文件共 119.03 MiB，首次获取四库合计约 35 秒。原件已就位时，导入和重读共 6.637 秒，进程峰值 RSS 661.92 MiB。4 个 H5AD 共 208.72 MiB，压缩 MTX 流式读取，未生成解压矩阵文件。每库写入临时 H5AD，通过重读后移到正式路径。完整环境与代码版本见 [依赖锁](python-packages.txt) 和 [资源记录](resource_profile.json)，逐库日志见 [import.log](import.log)。

从项目根目录重建。随机抽查种子为 20261121。

```bash
uv venv --python 3.12.13 .venv-m1
uv pip install --python .venv-m1/bin/python -r artifacts/m1/python-packages.txt
uv run --no-project --python .venv-m1/bin/python python -Bu scripts/import_gse245552.py \
  --manifest data/manifests/m0/sample_manifest_validated.csv --donor GSE245552:Pt18
```

[导入脚本](../../scripts/import_gse245552.py) 会复用本地一致的 GEO 原件，缺件按记录的 URL 获取。[M1 报表与日志 SHA](output_sha256.txt)

## M2 交接

其余 35 个 GSM 可从同一份 M0 样本表分批处理，复用逐库读取与 QC，保留 Pt08 的多个病灶。下一项验收是统一注释后，各谱系在 PT 与 LM 两侧同时有细胞支持的独立患者数。全队列统一 QC、程序评分、患者效应、模型比较和外部复现均待后续阶段。
