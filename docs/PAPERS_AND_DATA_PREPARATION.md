# 论文与数据准备清单

核验基线：2026-11-21。本清单与 `docs/CRLM_Research_Plan.md`、`docs/source_brief/说明.md` 配套；数据角色和冻结规则以研究方案为准。

**先准备 12 篇核心原始文献：8 篇服务于数据与统计设计，4 篇服务于首轮建模。** 候选机制、空间和实验文献按实际进入的分支追加。用户已经提供的 2020 年 Cancer Medicine PDF 保留即可，不用重复下载。`references/reading_manifest.csv` 保留完整 85 条参考记录，但这不是“需要立即下载 85 篇论文”的要求，更不是“85 篇全文已经下载”的声明。

参考目录包含 78 条非空 DOI 记录、77 个唯一 DOI；B35 与 M24 是同一篇 scPAFA 论文。部分 M 记录是官方代码或文档，本来就不需要 PDF。以下“本轮取得”描述的是此前研究过程中实际查阅的材料；“包内已有”只描述本项目 `starter/` 中确实携带的文件。二者不能互换。

## 1. 文件怎样准备

正式论文可保存为 `references/papers/B07_2024.pdf` 这类文件名；出版社或 PMC 的完整 HTML/XML 也可作为正式 Methods 的可检索来源，保存时保留标题、DOI、版本和源 URL。Supplementary Table、DOCX、XLSX、代码 ZIP 等放在 `references/supplements/<ref_id>/`，保留原始文件名。建议另记每个文件的来源 URL、取得日期、字节数和 SHA256；确认文件能打开，不能把下载失败的 HTML 页面改名当 PDF。

论文的 PDF、补充表、作者代码和数据元信息分别解决不同问题。PDF 用来核对假设、取材和统计方法；临床表与 sample metadata 用来决定真实患者配对；代码用于核对实现；表达矩阵与 feature/barcode 用来决定能否运行。**仅有 PDF 不足以重建可复现分析，下载了 PDF 也不等于作者代码或数据已取得。**

本项目已附事实型小证据、少量原始补充表、提取结果和审计脚本；未重新分发整批论文全文、完整表达矩阵或预训练权重。CSV 的 `full_paper_pdf_in_starter` 因此均为 false；B19 的两份表格 PDF 是原始补充表，不是论文全文。

## 2. 第一批：启动必需的 8 篇

“启动必需”指建立可信分析方案所必需的知识来源。离线文件检查和合成数据自检可先做；某篇仍缺完整 Methods 时必须登记缺口，不得用猜测填上原作者方法。PhenoAligner 原实现未取得会限制精确原算法复现，但无需让整个项目停在这里。

### 1）B07 — 发现集 GSE245552 的原论文

**Th17 Cells Secrete TWEAK to Trigger Epithelial–Mesenchymal Transition and Promote Colorectal Cancer Liver Metastasis**（2024，Cancer Research）。[DOI：10.1158/0008-5472.CAN-23-2123](https://doi.org/10.1158/0008-5472.CAN-23-2123)；[原始 GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245552)。

- **优先补齐：** 正式全文 PDF 或合法作者稿，尤其组织处理、富集/分选、QC、注释、治疗与取样时序的 Methods；若正文转引补充方法，也应一起取得。
- **本轮取得深度：** 仅摘要与官方数据记录；正式完整正文尚未取得。另已取得官方原始 S1/S2/S3，并完整核验 39 份公开 MTX 的数值与维度。这些文件审计不能代替论文全文阅读。
- **包内已有：** `starter/audit_data/papers/can-23-2123_table_s1_suppst1.docx`；`starter/audit_data/papers/GSE245552_supplement/` 中 S2/S3 原 DOCX；患者映射、feature 集合、逐矩阵历史审计 JSON 和真实 manifest。
- **必须解决的问题：** 公共 scRNA 是 16 患者、39 标本，13 对 PT–LM；不是按全研究 18 人或 130 份多组学材料计算。治疗、CD45 富集/分选与混合比例仍 unknown，不从 S2 的细胞产量推断。S2 的计数单位是 ×10⁶，不是 scRNA QC 后细胞数；S3 使用最后的 `gene` 列，不用带去重后缀的行名。

原表来源：[S1](https://doi.org/10.1158/0008-5472.25607328.v1)、[S2](https://doi.org/10.1158/0008-5472.25607325.v1)、[S3](https://doi.org/10.1158/0008-5472.25607322.v1)。已附原表无需重复准备。

### 2）B02 — PhenoAligner 与 GSE164522

**Immune phenotypic linkage between colorectal cancer and liver metastasis**（2022，Cancer Cell）。[DOI：10.1016/j.ccell.2022.02.013](https://doi.org/10.1016/j.ccell.2022.02.013)；[PubMed 原始记录](https://pubmed.ncbi.nlm.nih.gov/35303421/)；[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164522)。

- **优先补齐：** 正式全文 Methods、Supplemental information、个体临床/处理表，以及原 PhenoAligner 实现。重点核对参照细胞如何选取、KNN/距离定义、匹配阈值、检验统计单位和跨患者聚合。
- **本轮取得深度：** 摘要、出版社数据/代码声明、官方 GEO 与作者机构图解；正式完整 Methods 未取得。作者代码声明为 upon request，不能写成已经拿到原实现，也不能只凭二手图解复现后称“原算法”。
- **包内已有：** `starter/evidence/PhenoAligner_access_and_estimands.txt`、10 位公开患者身份快照、六组织完整列头匹配与表达尺度探针；未附完整论文或完整表达矩阵。
- **使用边界：** 公开四组织子集是 10 人/40 组织/163,347 细胞。公开表达是与 log1p(normalize10000) 数值兼容的处理后矩阵，不能逆变换后冒充 UMI。患者四组织交互与细胞表型对齐不是同一估计对象。

### 3）M24（同 B35）— 最直接的四组织方法前作

**Uncovering disease-related multicellular pathway modules on large-scale single-cell transcriptomes with scPAFA**（2024，Communications Biology）。[DOI：10.1038/s42003-024-07238-7](https://doi.org/10.1038/s42003-024-07238-7)；[PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11569158/)；[官方 PDF](https://www.nature.com/articles/s42003-024-07238-7.pdf)。

- **准备：** 正文、Supplementary information、[软件仓库](https://github.com/ZhuoliHuang/scPAFA)与[论文复现仓库](https://github.com/ZhuoliHuang/scPAFA_paper)；记录实际版本。特别读 CRC–CRLM 结果、Fig.4、四组织输入、程序定义和组间检验。
- **本轮取得深度：** 正式 PDF 全文，重点方法、CRC–CRLM 结果与官方复现入口；包内保留引用和研究判断，未附全文或作者仓库。
- **作用：** 决定创新边界和必要基线。它已经使用同一 GSE164522 的 10 人/40 组织/163,347 细胞；不能把它当新外部队列，也不能把复现旧程序称为首次发现。B35 与 M24 只下载一套。

### 4）B03 — 治疗异质性与 GSE178318 次要外部

**A single-cell atlas of liver metastases of colorectal cancer reveals reprogramming of the tumor microenvironment in response to preoperative chemotherapy**（2021，Cell Discovery）。[DOI：10.1038/s41421-021-00312-y](https://doi.org/10.1038/s41421-021-00312-y)；[PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8421363/)；[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178318)。

- **准备：** 正文 Methods、临床/治疗与样本映射补充表、QC/注释步骤、GEO 文件说明。实施时还需可连接到 barcode 的细胞注释；不能仅凭一份合并 MTX 推断身份。
- **本轮取得深度：** 已核关键结果、Methods 与图注；未声称逐字精读全部补充材料。包内已有原 GEO 元数据及队列审计，未附完整正文/矩阵。
- **作用与边界：** 6 位患者、6 对同期切除 PT–LM，3 位术前治疗、3 位未治疗；不是同患者治疗前后。小样本仅支持预设方向复现，不承担大规模外部分类性能结论。
- **版本核对：** 另有 [2023 年作者更正](https://doi.org/10.1038/s41421-023-00537-z)，更正作者 Rui-Min Xu 的姓名拼写；此次核验的更正没有提供新分析方法或新增患者。使用正式修订版本即可，此更正不另计为一篇必须精读的核心方法文献。

### 5）B46 — GSE203608 与人源共培养实验先例

**Mapping and modeling human colorectal carcinoma interactions with the tumor microenvironment**（2023，Nature Communications）。[DOI：10.1038/s41467-023-43746-6](https://doi.org/10.1038/s41467-023-43746-6)；[PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10689473/)；[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203608)；[作者处理对象记录](https://zenodo.org/records/7872684)。

- **准备：** Methods 的组织分选、样本标签、M-CSF 巨噬培养、类器官/CAF 共培养，Fig.1A/B，以及两版 cell annotation 的身份对应。需要表达时先选清一个版本。
- **本轮取得深度：** 核对官方全文关键方法/结果、原图临床与两份完整 metadata；未读取完整表达矩阵。包内已有 Fig.1、临床转录表、身份快照、配对支持量及去重审计。
- **作用与边界：** 真实体内子集 15 人/24 组织，只有 4 对 PT–LM、1 个完整四组织患者；该患者正常结肠的髓系只有 11 细胞。扩展 metadata 146,988 行已包含旧版 126,014 行全部，不得拼接。培养与技术库不计作新患者。

### 6）M09 — 先建立患者重复单位

**Confronting false discoveries in single-cell differential expression**（2021，Nature Communications）。[DOI：10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)。

准备完整 Methods 中的生物重复、伪重复与基线评估。此轮读取深度是原始摘要及关键方法/结果段落；正式实现前补读完整统计方法。它支持“以患者而非细胞作为独立重复”，不意味着任意小样本 pseudobulk 都自动校准。starter 有本项目透明患者聚合代码和审计材料，未附这篇全文。

### 7）M10 — 多患者、亚群状态的计数分支

**muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data**（2020，Nature Communications）。[DOI：10.1038/s41467-020-19894-4](https://doi.org/10.1038/s41467-020-19894-4)；[官方 vignette](https://bioconductor.org/packages/release/bioc/vignettes/muscat/inst/doc/analysis.html)。

准备完整 Methods、所安装版本的 vignette、sample × subpopulation 聚合与设计矩阵示例。本轮深入读取原文方法和官方 vignette，未把软件当作已经完成真实数据重分析。计数分支仅用于可追溯整数 counts；GSE164522 的处理后 log 表达走独立的程序分数分支。包内未附原论文或 muscat 全软件环境。

### 8）M13 — 固定程序的 rank 评分

**UCell: Robust and scalable single-cell gene signature scoring**（2021，Computational and Structural Biotechnology Journal）。[DOI：10.1016/j.csbj.2021.06.043](https://doi.org/10.1016/j.csbj.2021.06.043)。

准备原始评分 Methods 和安装版本的官方实现，核对 `maxRank`、负向基因、未测基因和排序背景。本轮已读原始方法记录及官方实现，未逐页声称全文精读。跨数据评分要同时冻结可测候选集合与排序的背景基因宇宙；只有 signature 基因取交集仍不够。包内有覆盖审计与候选种子，未包含 UCell 论文全文。

## 3. 第一批的后 4 篇：首轮建模前读

这 4 篇使首批总量为 12。先实现患者分组、简单基线与不泄漏的验证，再决定 scGPT 是否值得保留。

| 顺序 | 文献与原始入口 | 应准备的材料、用途及实际阅读深度 |
|---|---|---|
| 9 / M20 | **Attention-based Deep Multiple Instance Learning**，2018，ICML/PMLR；[官方正文及 PDF 入口](https://proceedings.mlr.press/v80/ilse18a.html)。此条无 DOI。 | 原始 attention pooling 公式、bag/instance 定义、损失；本轮核对原始架构和公式。starter 的透明 MIL 是项目原型，不是作者仓库的完整复现；attention 不等于因果解释。 |
| 10 / M02 | **Incorporating hierarchical information into multiple instance learning for patient phenotype prediction with single-cell RNA-sequencing data**，2025，Bioinformatics；[DOI 10.1093/bioinformatics/btaf241](https://doi.org/10.1093/bioinformatics/btaf241)；[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12261414/)；[官方代码](https://github.com/minhchaudo/hier-mil)。 | 正式主文、方法/实验、层级划分、PCA/scGPT 输入与官方代码版本。本轮读作者机构正式 PDF 和代码。它是“细胞类型层级＋scGPT＋MIL”的直接先例，组合本身不能作为创新。 |
| 11 / M07 | **scGPT: toward building a foundation model for single-cell multi-omics using generative AI**，2024，Nature Methods；[DOI 10.1038/s41592-024-02201-0](https://doi.org/10.1038/s41592-024-02201-0)；[官方仓库](https://github.com/bowang-lab/scGPT)。 | 主文预训练 Methods、M19/M25–M27 对应代码、实际 checkpoint/vocab/config/许可。本轮读原文/预训练方法并审计固定 commit `cebd6fae655b9c585a4807daa3ac31bb764f06b4` 的相关实现。starter 没有模型权重；冻结不能自动排除预训练重叠。 |
| 12 / M08 | **Zero-shot evaluation reveals limitations of single-cell foundation models**，2025，Genome Biology；[DOI 10.1186/s13059-025-03574-x](https://doi.org/10.1186/s13059-025-03574-x)；[官方复现仓库](https://github.com/microsoft/zero-shot-scfoundation)。 | 原始评估设计、PCA 等基线、输入与任务边界。本轮核读方法/评估与官方复现入口；不必为本项目下载其全部 benchmark。它支持公平对照，不支持预判本项目 scGPT 必然失败。 |

以上四篇原文/作者完整仓库未装入 starter；已包含引用、实现审查结论和本项目自有原型。不要把原型自检通过写成“复现作者论文结果”。

计数差异分析实际进入时再读 S04 edgeR v4（[10.1093/nar/gkaf018](https://doi.org/10.1093/nar/gkaf018)）及所安装版本指南；组成分析、cNMF、替代 MIL 等只在预先选定该分支时读对应条目。`reading_manifest.csv` 中大量 `02_first_model_conditional` 是按需方法储备，不是要求把所有方法都做一遍。

## 4. 候选程序、空间与实验：进入相应分支前追加

这里的“后读”是相对前述设计核对而言。**准备冻结 R/C/A 候选程序时，应先补读 B06、B19；准备实验时再补实验对应文献。** 不应等外部结果显著后再挑选能支持它的论文。

| 优先级/条目 | 原始论文或来源 | 要取得什么；包内已有与关键限制 |
|---|---|---|
| 候选冻结前 / B06 | **Colorectal Cancer Metastases in the Liver Establish Immunosuppressive Spatial Networking between Tumor-Associated SPP1+ Macrophages and Fibroblasts**，2023；[DOI 10.1158/1078-0432.CCR-22-2041](https://doi.org/10.1158/1078-0432.CCR-22-2041)；[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9811165/)。 | 已核关键正文/方法，取得官方 S1–S8 XLSX；原 XLSX 在 `starter/literature_seeds/literature_program_sources/B06_supplement/`，无需重下。读 S1 抗体限制、S4 DE、S5 富集交集、S6 非癌背景。S5 是报告交集，不是完整 Reactome 通路；该研究无配对 PT。 |
| 候选冻结前 / B19 | **CXCL9:SPP1 macrophage polarity identifies a network of cellular programs that control human cancers**，2023；[DOI 10.1126/science.ade2292](https://doi.org/10.1126/science.ade2292)；[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10755760/)。 | 已核正文关键段，取得并视觉核验原 Table S4/S5。两份原表 PDF、Mac_SPP1/Mac_CXCL10 完整块、561 条巨噬相关记录及提取脚本均在 `starter/literature_seeds/`。仍需正文/Methods；HNSCC 程序不是 CRLM 特异签名，相关统计不是因果权重。 |
| 使用 225 前 / B05 | **Single-cell and spatial transcriptome analysis reveals the cellular heterogeneity of liver metastatic colorectal cancer**，2023；[DOI 10.1126/sciadv.adf5464](https://doi.org/10.1126/sciadv.adf5464)；[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10275599/)。 | 本轮取得原补充 PDF 和 code.zip；包内只有 S2/S3/S4 转录与证据、身份核查。必须另取得代码引用的 **sampleinfo.csv** 或等价作者身份映射并解释 s0816。该文件本轮未取得；不凭日期猜别名，不删除 0816 使分母看起来一致。 |
| 正常器官解释 / B18、B21 | B18 **Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches**，2022，[10.1016/j.cell.2021.12.018](https://doi.org/10.1016/j.cell.2021.12.018)；B21 **Cross-tissue immune cell analysis reveals tissue-specific features in humans**，2022，[10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)。 | 读肝驻留身份及真实多组织供者 Methods。本轮核关键正文/方法；包内是正常图谱 metadata 支持审计，无完整表达。不能把全图谱总供者当肝–结肠配对 N。 |
| 功能实验前 / B08、B11、B38 | B08 **SPP1+ macrophages cause exhaustion of tumor-specific T cells in liver metastases**，2025，[10.1038/s41467-025-59529-0](https://doi.org/10.1038/s41467-025-59529-0)；B11 **Liver metastasis restrains immunotherapy efficacy via macrophage-mediated T cell elimination**，2021，[10.1038/s41591-020-1131-x](https://doi.org/10.1038/s41591-020-1131-x)；B38 **SPP1high macrophage-induced T-cell stress promotes colon cancer liver metastasis through SPP1/CD44/PI3K/AKT signaling**，2025，[10.1136/jitc-2025-012330](https://doi.org/10.1136/jitc-2025-012330)。 | 准备具体细胞来源、共培养/迁移、干预与恢复、抗体、剂量时程、阴性对照的 Methods/补表；本轮核关键机制与方法，不等于全部实验条件已经可直接照做。B38 四组织人数据仅 reasonable request，未落实为现成公开验证集。 |
| 空间分支 / B27 | **Evidence of endothelial anergy and T-cell exclusion in microsatellite stable colorectal liver metastases revealed through human hemi-liver perfusion and spatial transcriptomics**，2026；[10.1186/s12943-026-02732-4](https://doi.org/10.1186/s12943-026-02732-4)。 | 本轮核 Methods/结果/数据声明；4 位 MSS、replacement 型、术前治疗 CosMx 患者不能代表所有 CRLM。GSE312757 大对象仅入口/大小核查，未读约 7.68 GB 表达对象；先核面板、部位与可测基因再下载。 |
| 最新机制边界 / B44、B45、B47 | B44 **Senescent-Like Neutrophils Shape Angiogenic Immunosuppressive Niches in Colorectal Cancer Liver Metastasis**，[10.1158/2159-8290.CD-25-0820](https://doi.org/10.1158/2159-8290.CD-25-0820)；B45 **APOC1⁺IFI30⁺ TAMs orchestrate colorectal cancer liver metastasis by coupling cholesterol biosynthesis with secretion of immunosuppressive factors**，[10.1007/s00018-026-06338-3](https://doi.org/10.1007/s00018-026-06338-3)；B47 **Macrophage PABPC4-SPP1 Axis Orchestrates Immunosuppression in Colorectal Cancer**，[10.1002/advs.77361](https://doi.org/10.1002/advs.77361)；均为 2026。 | 已核官方全文关键 Methods/结果及相关补表，逐篇边界见 CSV 和研究报告。B44 人样本中性粒富集重混，不能复现天然组成；B45 复用旧队列，不能当新 46 人验证；B47 含肝转移阴性表述，需保留不支持部分。包内是研究备忘，未附整批全文/补充。 |
| 互作分支 / S01–S03、M23 | [S01 资源/方法比较](https://doi.org/10.1038/s41467-022-30755-0)、[S02 LIANA+](https://doi.org/10.1038/s41556-024-01469-w)、[S03 NicheNet](https://doi.org/10.1038/s41592-019-0667-5)、[M23 MultiNicheNet](https://doi.org/10.1101/2023.06.13.544751)。 | 读原始 Methods、配体受体资源版本、患者重复处理与输出含义。M23 当前参考表标为预印本。推断分数不能证明分泌、受体结合或下游功能；进入此分支再准备，不作为第一周必须装齐的软件。 |

B06/B19 原表的公开来源及哈希在 `starter/literature_seeds/literature_program_seeds_README.txt` 和 `literature_program_sources/seed_manifest.json`。种子共 310 条“集合—基因”记录、293 个原始符号、10 组，只是可追溯候选材料；不等于已经冻结或验证的本项目 signature。

## 5. 数据准备比新增论文下载更早的检查

首轮只准备发现集，冻结后按计划打开外部分析。下表是将来实际运行需要的文件，不表示这些大数据已经放入项目。

| 队列与角色 | 必须实际取得并核对 | 已有证据与尚缺材料 |
|---|---|---|
| **GSE245552：唯一发现/开发集** | 39 GSM 的 MTX、features、barcodes，与 `starter/audit_data/GSE245552_sample_manifest_audited.csv` 一一对应；S1/S2/S3；预处理说明 | 历史 39 MTX 全文件审计已完成；压缩 MTX 合计约 1.06 GB，但矩阵本身不在包内。33 肿瘤为 29,569 genes，6 邻近组织为 29,922；缺失 353 genes 含所有 13 个 MT genes。肿瘤 pct_mt 必须 NA，不能从删过的特征表算出 0 后当作高质量。 |
| **GSE164522：主要外部四组织确认** | 官方六组织表达、作者身份/注释；只选择预定 PT/LM/NC/NL；保留原处理尺度与真实可测宇宙 | 包内六表头/身份匹配已核，共 257,843 列；四组织 10 人/163,347 细胞。完整表达未附。Myeloid=Mono/Mac/DC 时，50 细胞/组织门槛仅 5 个完整四组织患者。外部 metadata 已被本轮查看，需保留暴露日志。 |
| **GSE178318：次要方向复现** | 合并表达及 feature/barcode、患者/组织/治疗映射、可验证细胞注释 | 6 对 PT–LM；本轮输入 140,281 barcode 与论文 QC 后 111,292 不同，不能混用。完整表达未附；保留 GEO 原数据使用说明。 |
| **GSE203608：次要独立来源候选** | 选定一个版本的表达与匹配 metadata；明确只纳入 in vivo；处理技术标签/培养样本 | 身份 metadata、4 对的临床图和髓系支持已核；完整矩阵数值未核。不得把旧/新 metadata 拼接，也不把 53 GSM 当 53 患者。 |
| **GSE225857：暂缓配对主分析** | sampleinfo.csv 或等价身份说明；正文与 S2/S3/S4/代码对应；RNA 与空间病例分开 | 6 临床代码与公开 7 日期代码不一致，s0816 未解。BD Rhapsody 是全转录组加 7 AbSeq，并非 7 个 RNA 基因；空间 4 PT/2 LM 不配对，也不与 scRNA 同人。 |
| **GSE298084：小型探索；GSE335552：空间可测性** | 前者真实 3 人/8 文库及围瘤位置；后者每切片 panel、segmentation 和样本身份 | 前者不是 18 人 scRNA；后者 3 患者/4 Xenium 切片，无 PT。包内已有实际面板：4 切片共同 309 基因；未测候选不能补 0 当阴性。不是首轮必须下载的数据。 |

具体官方 URL、许可/访问限制、收集来源与重叠风险见 `starter/evidence/data_sources_D01_D10.txt`（文件名沿用，正文已含 D01–D11）、`dataset_audit.csv`、`cohort_provenance_audit.csv`。不同 accession 只能支持独立来源审查，不能绝对保证没有跨研究重复患者。

无须现在下载所有正常图谱表达、FASTQ、全部空间图像、所有 foundation-model 权重。先完成身份、测量覆盖、版本与资源检查；受控数据当前不可取得时，维持不可用状态，不加入已落实样本量。

## 6. 交给无联网 Agent 时

1. 把**实际项目文件**放进目标工作目录，包括 `docs/source_brief/说明.md`、研究方案、这份清单、`references/reading_manifest.csv` 和完整 `starter/`。用户原 2020 PDF 已放在 `references/papers/B01_2020_classifier.pdf`；请 Agent 直接核对并读取；不要只给一个聊天中的文件名。
2. 把已合法取得的首批论文/Methods、原补充表、必要作者代码一并带入。请 Agent 先列出真实存在的路径、打开文件核对标题/DOI/页数或内容，再更新阅读状态。以目标环境可见文件为准，不依赖另一个聊天或工作区可能保存的上下文。
3. 优先运行 starter 的默认离线证据复核和合成示例。它们可以检查所附小证据与代码流程，不能重新证明未附的 1 GB 表达矩阵正在本机，也不能产生真实生物学结论。
4. 缺 B07 全文、B02 原 Methods/实现、225 身份文件、scGPT 权重等，逐项记 missing。仅跳过依赖该材料的工作，保留可完成的审计、基础框架和统计设计；不得用推测替代临床或作者实现。
5. 真正跑数据前，由用户准备指定真实矩阵、对应元数据和锁定软件环境。只具备论文 PDF 的离线环境不足以运行真实研究。

## 7. 交给可联网 Agent 时

1. 同样先检查本地文件；已有的 S1/S2/S3、B06 XLSX、B19 表格 PDF 和原 2020 PDF无需重复下载。先读 `reading_manifest.csv` 的 `first_batch`、`required_material`、`missing_or_blocking`。
2. 用 DOI、出版社、PMC/Europe PMC、作者官方仓库逐一取得缺件。保留源 URL、版本/commit、取得日期、大小与 SHA256。URL 指向摘要/元信息时，要继续明确它不是正式全文；403、登录或受控数据限制只记录，不能假装已取得。
3. 先拿小文件：临床表、sample metadata、feature/barcode、面板、代码配置。确认患者身份与测量约束后，才下载当前阶段真正需要的表达文件；不自动批量下载所有参考文献的数据。
4. 模型权重要单独核对 checkpoint、vocab、预处理、许可与预训练来源。当前仓库源码不能独自证明某个发布权重的全部训练病例；把未排除的预训练重叠保留为限制。
5. 每篇完成后只更新实际阅读的深度。例如“取得 PDF，读 Methods/S1”与“逐图核查、复跑作者结果”是不同状态。没有作者回复时，不把 upon request 的数据/代码写成已取得，也不自行发邮件。

## 8. 完成准备的验收条件

- 首批 12 篇逐项有准确标题、年份、DOI/官方入口及实际文件状态；B01 原附件保留，B35/M24 去重。
- B07 全文、B02 原 Methods/实现、225 身份映射等缺口明确可见；已有引用或审计摘要没有被标成“全文已在包内”。
- 所有 starter 路径在新项目中能打开；数据角色保持 245 发现、164 主要外部、178/203608 次要外部、225 映射未解时暂缓。
- 程序种子、作者 marker、富集交集与因果机制证据分开标识；冻结候选前读相关原文，不在外部结果揭晓后选择证据。
- 第一个真实数据分析运行能追溯到实际文件、患者 manifest、输入尺度、测量基因、预处理版本、患者划分及随机种子；这比“论文下载数量”更能说明准备是否完成。

CSV 的字段 `starter_relative_paths` 是相对于本项目根目录的 `starter/...` 路径；`proposed_local_paper_path` 只是推荐落盘位置，**不是当前已存在文件的声明**。具体逐条阅读深度、原文链接、材料要求和未决限制以该 CSV 为可更新记录。
