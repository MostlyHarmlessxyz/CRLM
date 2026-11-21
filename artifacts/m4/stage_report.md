# M4 · 患者分折的组织来源分类

在 GSE245552 的五折留患者预测中，固定训练折 PCA 后取细胞均值，T/NK、髓系和 B/Plasma 的患者组织级 AUROC 分别为 1.000、0.988、0.867。髓系的线性 attention 和 MixMIL 均为 0.994；MixMIL 相对于同表示均值头的差值为 +0.006，按共同患者重抽样的 95% 区间为 0.000 至 +0.031。复杂模型没有表现出跨谱系一致的增益。[完整指标](patient_metrics.csv) · [直接前作比较](direct_prior/stage_report.md)

| 固定主种子下的 AUROC | T/NK | Mono/Mac/DC | B/Plasma |
|---|---:|---:|---:|
| 程序 Elastic Net | 0.513 | 0.828 | 0.175 |
| 组成 Elastic Net | 0.492 | 0.864 | 0.467 |
| PCA32 均值 | 1.000 | 0.988 | 0.867 |
| PCA32 线性 attention | 0.862 | 0.994 | 0.683 |
| scGPT512 均值 | 0.923 | 0.941 | 0.808 |
| scGPT512 线性 attention | 0.667 | 0.959 | 0.442 |

稀疏处理后有 33 个合格肿瘤文库、89 个标本×谱系袋、18,610 个固定抽样细胞。各谱系可贡献分类指标的患者分别为 16、15、15，其中有完整 PT–LM 对的为 12、11、8。折按患者分配，PCA、特征和缩放仅在训练患者上拟合，同部位多病灶先等权合并。模型识别已经取到的组织。训练和统计使用 CPU，冻结编码使用本机 Intel Arc A750。[分折](split_manifest/eligible_samples_by_fold.csv) · [抽取细胞](selected_cells.csv.gz) · [预定比较](benchmark_design.md)

共登记 390 次拟合并得到 2,279 条折外袋预测。269 次收敛，120 次完成预定步数，T/NK 的一个 scGPT attention 次要种子在第 0 折达到 500 次迭代上限。该次保留失败状态；主种子结果完整。三谱系的 15 个迁移头均在整个开发集重拟合并收敛。[各次运行](benchmark_runs.csv) · [折外预测](oof_predictions.csv) · [全开发拟合](full_data_models.csv)

贡献定义为每个细胞对袋预测 logit 的带符号加项，attention 为池化权重。逐袋移除绝对贡献最高的一组细胞，同时用亚型和数量匹配的随机移除作为参照。可配对的 350 组折外跨种子比较中，PCA 和 scGPT 的带符号贡献秩相关中位数分别为 1.000 和 0.594，绝对贡献前 10% 集合的 Jaccard 中位数分别为 1.00 和 0.25。scGPT 的细胞归因对种子敏感。[稳定性](attribution_stability.csv) · [PCA 移除试验](PCA/attribution_stress.csv) · [scGPT 移除试验](scGPT/attribution_stress.csv)

12 个真实细胞中，A750 与相同权重 CPU 结果的最小余弦相似度为 0.99999976，小批测试峰值设备分配为 656 MB。33 个标本缓存总共 35.8 MB，各文件读取、细胞身份、向量宽度、有限值和单位范数均已核对。11,008 个细胞超过 1,199 个表达 token 的预算，采用冻结的逐细胞抽样。开发库是否曾进入预训练语料尚无逐患者清单可查。[权重及来源](resource_profiles/scgpt_sources.md) · [小批实测](resource_profiles/scgpt_profile.json) · [缓存记录](resource_profiles/scgpt_cache_manifest.csv)

复算顺序为 `benchmark_discovery.py --representation PCA`、`encode_scgpt.py --device xpu --batch 4`、`benchmark_discovery.py --representation scGPT`，然后按直接前作报告运行三个 `benchmark_prior_mil.py` 分支，最后运行 `fit_transfer_models.py` 和 `summarize_benchmarks.py`。入口都在 `scripts/`，运行环境分别见 [CPU 环境](python-packages.txt)、[XPU 环境](resource_profiles/xpu-packages.txt) 和 [直接前作环境](direct_prior/extra_packages.txt)。各模型输入来自 M2 的真实表达和 M3 的预定义文献程序。
