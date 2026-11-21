GSE245552 原始补充资料 S2/S3（D07k）

S2: DOI10.1158/0008-5472.25607325.v1
https://ndownloader.figshare.com/files/45666639
S3: DOI10.1158/0008-5472.25607322.v1
https://ndownloader.figshare.com/files/45666636
原文件大小、SHA256、对应API URL见同目录GSE245552_S2_S3_provenance.json。
原DOCX存于../papers/GSE245552_supplement/，转录CSV在本目录亦有副本；
离线脚本读取audit_data主目录的同内容CSV。

S2共有39标本、16患者，结构为16PT/17LM/4相邻结直肠/2相邻肝。
CD45+/CD45−数值原单位为×10^6，已从DOCX XML核实两处superscript6。
这些百万级数值不是scRNA条形码数、QC细胞数或任何谱系有效N；不能推断
细胞计数所在实验阶段或FACS富集/混合比例，也不把S2的PT1/LM1编号硬配GEO_T序号。
S2写Colon tumor不能覆盖S1的精确原发部位，例如其中包含临床直肠病例。

S3五表共600marker行、390唯一gene、30个作者群；390全部在公共肿瘤feature中。
必须使用最后的gene列；首列是R行名，含PLVAP1/HSPG22/PECAM17等去重后缀。
转录CSV明确标识作者cell-level p值/FC，不能当本项目患者效应或独立复现。
原作者cluster名及marker只作为注释先验，不能直接认定跨队列同名状态完全相同。

离线重跑（在CRLM_starter目录）：
python audit_data/audit_GSE245552_supplement_tables.py
默认核心离线审计reproduce_core_audit.py亦读取这些原件检查单位、身份和gene覆盖。
