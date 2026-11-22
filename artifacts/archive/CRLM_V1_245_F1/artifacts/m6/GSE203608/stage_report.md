# UPenn 独立次要来源

在预登记的 GSE203608 体内组织中，髓系 Mac_SPP1 的四个患者差为 +0.00795、+0.01228、−0.00561、−0.03458，均值 −0.00499；B/Plasma UPR 的两位合格患者分别为 −0.00984 和 +0.03331。两项的双侧精确 sign test 原始 p 均为 1，均未满足预期方向。去除原作者图示 Stage III 而注释含 LiverMet 的 49 号后，髓系均值为 +0.00487，2/3 同向。[逐患者数据](patient_effects/patient_contrasts.csv) · [临床敏感性](patient_effects/exclude_patient49.csv)

输入为原版 126,014 个细胞的计数矩阵、原版细胞及基因注释，逐行对应历史身份表。读取完全部 315,102,049 个整数坐标，33,694 个作者特征中 33,626 个有唯一基因符号。通过作者 QC 的体内标本含 59,661 个细胞、24 个标本、15 个患者代码，四对解剖 PT–LM。Mac_SPP1 的 160 个冻结成员有 149 个可测，UPR 为 106/107。门槛 20、50、100 个细胞时，髓系始终四对，B/Plasma 分别为三、二、一对。[患者与标本](coverage/sample_manifest.csv) · [基因覆盖](coverage/measurement_coverage.csv) · [原始研究](https://doi.org/10.1038/s41467-023-43746-6)

在 M4 CPU 环境运行 `scripts/fetch_gse203608.py` 取得原始表达，再运行 `scripts/validate_gse203608.py`。矩阵和注释来自同一原始版本，扩展培养样本不进入分析。真实分析用时 38 秒，峰值 RSS 1.05 GB。三个输入的 SHA 及坐标规模可用于复算。[资源与来源](resource_profile.json) · [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203608)
