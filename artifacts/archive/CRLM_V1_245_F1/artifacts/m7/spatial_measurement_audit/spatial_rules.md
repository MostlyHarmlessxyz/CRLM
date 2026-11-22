# GSE335552 实测空间读出

四张 Xenium 切片来自三位 CRLM 患者，M131 提供肝转移及邻近肝，其余两位仅有肝转移。逐切片将原始细胞×特征 H5 与细胞质心按完整条码一一匹配，保留 Gene Expression 探针。几何单位为微米。[GEO 数据与八个输入附件](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335552) · [10x 坐标定义](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs)

对两个冻结程序列出各切片实测成员。代表性标记的主要问题为，同一局部区域内 SPP1 阳性与阴性髓系细胞到最近 CD3/CD8 细胞的距离是否有差异。髓系须在 CD68、CD163、C1QA、CD14、LYZ 中测到至少两个，且含前三个之一。CD8 T 的定义为 CD3D/E 与 CD8A/B 各至少一个。两者均要求 EPCAM 未检出，同时满足两种谱系规则者单列。SPP1 阳性阈值为至少 1 个转录本。

按 200 µm 方格比较，每个纳入格中两组髓系各至少两个。格内均值差为阳性减阴性，跨格按两组较小的细胞数加权。负值表示 SPP1 阳性髓系与 CD3/CD8 细胞间距较小。199 次格内标签置换保留局部密度与阳性数，给出置换参照范围。另按 100/400 µm 方格和至少 2 个 SPP1 转录本复算。每张切片分别报告，独立患者按 3 人计。区域无组织学分区时只解释面板上测到的局部关系。
