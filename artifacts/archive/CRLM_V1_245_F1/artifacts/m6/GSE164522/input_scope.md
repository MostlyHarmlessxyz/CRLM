# GSE164522 输入与身份

GEO 的四个公开表达文件分别为 PT（原发灶）、MT（肝转移灶）、PN（邻近结肠）及 MN（邻近肝）。每份按基因行、细胞列读取，保留作者非负连续值。V1 的四组织分析范围是该来源的 10 位公开患者（patient08–patient17）。原文报告的其他病例、PBMC、淋巴结和复用这些病例的后续论文不进入这一独立 N。[原研究](../../../references/papers/B02_2022.pdf) · [GEO 元数据](../../../starter/audit_data/geo/GSE164522_metadata_only.soft)

作者完整 metadata 有 257,843 个不同条形码，其文件 SHA 与既有原始记录完全相同。作者 `patient` 确定供者，`tissue` 确定组织，`sample` 确定供者内的取材标本。表达列使用 R 写出的点号形式，metadata 条形码按已发表交付格式把连字符映射为点号后与各组织匹配，两个版本同时保存。四组织共 163,347 个细胞、42 个作者标本。patient09 的肝转移标本 1、2 各有 645、282 个细胞，邻近肝标本 1、2 各有 4,492、2,687 个细胞；合格的编号标本在同一患者同一组织内等权。B02 对该队列报告入组时未治疗，逐例采集日期仍未知。[作者标本索引](coverage/sample_manifest.csv) · [细胞与标本映射](coverage/cell_identity.csv.gz)

作者字段 `celltype_major` 用于外部宽谱系映射。髓系只合入 Monocyte、Macrophage、DC，MAST 与 ILC3 各自保留；Bcell 包含 B 和浆细胞。报告同时保存 `celltype_global`、`celltype_sub`，便于按冻结规则检查主要类别与上层注释不一致的细胞和带跨谱系标志的髓系组。B02 正文记载 CRLM 细胞经 CD45 阳性分选。组成分母包括该标本进入作者矩阵的全部细胞，也保留其作者标注的 CD45− 细胞。GSE245552 使用 CD45 阳性与阴性细胞重混，两种方案的比例仅在各自捕获范围内解释。[标注映射](coverage/author_lineage_mapping.csv) · [B02 正文](../../../references/papers/B02_2022.pdf)

GSE164522 与 GSE245552 来自不同研究团队和临床来源。跨研究审计记录了 164 与其后续复用研究、245 与同队列组学的关系；公开作者代码并不提供逐人跨医院链接。既有审计已查看 164 的部分已发表亚组比例。正式外部检验随冻结记录标明该暴露。[独立来源审计](../../../starter/evidence/cohort_provenance_audit.csv) · [本阶段访问记录](external_access_log.csv)
