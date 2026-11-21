# GSE245552 身份与输入

本次在运行机将 [样本表](../../data/manifests/m0/sample_manifest_validated.csv) 逐 GSM 对到作者 GEO SOFT、逐库审计 JSON 和带原 URL 的清单。39 个 GSM 的标本号、来源组织、患者、输入 barcode 数和三项文件 URL 对齐。患者×组织的取材数与原 S2 一致，年龄、性别、原发位置、MSS 与 S1 一致。

| PT | LM | NC | NL | 患者 | PT–LM 完整患者 | 四组织完整患者 |
|---:|---:|---:|---:|---:|---:|---:|
| 16 | 17 | 4 | 2 | 16 | 13 | 1 |

13 对来自 Pt01、Pt02、Pt03、Pt04、Pt05、Pt06、Pt07、Pt08、Pt10、Pt11、Pt15、Pt17、Pt18。Pt09、Pt13、Pt16 缺一侧肿瘤组织。Pt08 有两份 PT、两份 LM，在配对推断中占一人；Pt16_T1 的 GEO 来源为肝转移。S1 的 Pt19、Pt20 无对应 GSM。各谱系、各程序最终可用患者数取决于后续覆盖。

Pt18 四组织为 GSM7844847/PT、GSM7844848/LM、GSM7844849/NC、GSM7844850/NL。共 22,179 个输入 barcode，一个 donor_uid。病灶/标本保留作者标本号，sample_uid 使用 GSM，分折与患者重抽样使用 donor_uid。

## 治疗、捕获与测量口径

B07 正文第 1353–1354 页报告该队列为同步肝转移、MSS，患者此前未接受抗肿瘤治疗。作者报告 CD45 阳性与阴性经磁选后按 2∶1 重混。S1/GEO 的逐人治疗字段仍为空，样本表沿用其来源标记 unknown；分析引用队列治疗史时单列 B07 依据。实际回收比例及各组织逐例取材日未知。组成估计限定于捕获细胞。解离条件存在来源差异，正文为胶原酶 II、20 分钟，GEO 记录试剂盒、30 分钟。

S2 产量单位为百万细胞。S3 原 DOCX 五表的 600 行组别与末列 gene 已对回 marker CSV，共 390 个原始符号。此表用于注释参考。历史 39 库共 190,190 个输入 barcode。肿瘤交付轴 29,569 个基因，邻近轴 29,922，肿瘤缺全部 13 个 MT 基因。历史 raw_counts 声明指作者交付的已筛选整数计数；本次 M0 复核来源及历史审计，实际选中库的稀疏数据由 M1 核验。

[临床与来源材料索引](input_inventory.csv) · [独立患者表](../../data/manifests/m0/donor_manifest_validated.csv) · [validator 结果](../../data/manifests/m0/manifest_validation.json)
