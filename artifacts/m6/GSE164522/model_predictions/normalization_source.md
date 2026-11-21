# GSE164522 肿瘤分类输入

Liu 等在 STAR Methods 的 “Single-cell RNA-Seq data processing” 中说明，10x 液滴数据经细胞质控后按文库总量归一到每细胞 10,000，再调用 `scanpy.pp.log1p()`。文中的回归和 PCA 在此之后进行。GEO 样本 CRLM08–17 的平台为 10x 5′。公开肿瘤表达表的 74,888 个细胞键均可与原 metadata 唯一匹配，均为 16 bp 碱基条码加文库后缀。GEO 的通用说明将补充表称为 UMI counts，实际交付为非负连续数值。[论文，STAR Methods e3](https://doi.org/10.1016/j.ccell.2022.02.013) · [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164522)

PCA 采用上述作者方法和冻结的逐标本数值判据。对每个非空细胞在作者完整特征轴求 `sum(expm1(X))`，其中位数须在 9,900–10,100，且至少 95% 的细胞位于 9,500–10,500。通过后直接使用作者表达层，按符号匹配训练期 512 个 HVG，至少 95% 可测才投影，未交付的列取训练中心。PCA 记录到的实测零、缺列和中心占位位置数分别列表。[M5 规则](../../../m5/validation_rules.md) · [PT 尺度](../coverage/PT_scale_check.csv) · [LM 尺度](../coverage/MT_scale_check.csv)

程序模型使用该来源四组织共同可测的基因轴和冻结的正向基因集，分母采用达到冻结覆盖规则的可测成员。scGPT 在作者交付值的正表达基因中按冻结实现分箱、截断并编码，保存每袋词表覆盖与截断细胞数。[程序覆盖](../coverage/candidate_measurement.csv) · [scGPT 开发配置](../../../m4/resource_profiles/scgpt_sources.md)
