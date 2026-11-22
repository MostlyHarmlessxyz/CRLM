# 冻结分类器在 GSE164522 的组织来源识别

T/NK 的 scGPT 均值池化在 10 位患者中的 AUROC 为 0.890，10 对患者均给出高于原发灶的肝转移灶 logit。髓系上，PCA 均值池化与注意力池化的 AUROC 均为 0.792。B/Plasma 的固定程序模型 AUROC 为 0.230，10 对患者中有 2 对的肝转移灶 logit 较高。三个谱系的 scGPT 注意力相对均值池化的 AUROC 增量分别为 −0.010、−0.069 和 −0.180，其配对区间均覆盖零。[患者指标](patient_metrics.csv) · [配对模型差](paired_method_contrasts.csv)

## 患者层结果

使用冻结版本 `72ee86b`。患者同一组织的各合格病灶按 logit 等权平均。区间来自 10,000 次按患者重抽样，阈值固定为 LM 概率 0.5。表中 AUROC、AP、balanced accuracy 使用该方法有输出的患者×组织，最后一列仅使用 PT 和 LM 都有结果的患者。髓系覆盖 9 位患者、17 个患者×组织、8 个完整患者对；T/NK 与 B/Plasma 各覆盖 10 位患者、20 个患者×组织、10 个完整患者对。[患者×组织 logits](patient_site_logits.csv) · [患者内差值](patient_paired_logits.csv)

| 谱系 | 冻结分类器 | AUROC（患者 bootstrap 95% 区间） | AP | Balanced accuracy | LM logit 更高的配对 |
|---|---|---:|---:|---:|---:|
| T/NK | 程序 Elastic Net | 0.590（0.340–0.820） | 0.684 | 0.500 | 6/10 |
| T/NK | scGPT mean | 0.890（0.790–1.000） | 0.918 | 0.700 | 10/10 |
| T/NK | scGPT attention | 0.880（0.680–1.000） | 0.829 | 0.700 | 9/10 |
| Mono/Mac/DC | 程序 Elastic Net | 0.611（0.333–0.840） | 0.736 | 0.521 | 5/8 |
| Mono/Mac/DC | PCA mean | 0.792（0.508–1.000） | 0.860 | 0.715 | 6/8 |
| Mono/Mac/DC | PCA attention | 0.792（0.569–0.963） | 0.837 | 0.688 | 7/8 |
| Mono/Mac/DC | scGPT mean | 0.639（0.397–0.889） | 0.756 | 0.549 | 6/8 |
| Mono/Mac/DC | scGPT attention | 0.569（0.333–0.827） | 0.711 | 0.528 | 7/8 |
| B/Plasma | 程序 Elastic Net | 0.230（0.060–0.440） | 0.398 | 0.500 | 2/10 |
| B/Plasma | scGPT mean | 0.530（0.300–0.760） | 0.646 | 0.550 | 5/10 |
| B/Plasma | scGPT attention | 0.350（0.120–0.590） | 0.457 | 0.500 | 3/10 |

## 覆盖与拒判

两个肿瘤组织共 74,888 个作者保留的 CD45⁺ 细胞、21 个标本及 10 位患者。每个标本×谱系至少需要 50 个细胞，固定抽样上限 256。三类谱系产生 58 个合格袋、13,551 个唯一细胞。五种方法始终使用同一袋的同一批细胞。21 个标本 × 3 个谱系 × 5 种方法共 315 次判定；208 次产生预测、25 次因细胞不足拒判、82 次因训练 HVG 覆盖不足拒判。315 行均保留。[完整袋记录](bag_predictions.csv) · [细胞名单](selected_cells.csv.gz) · [覆盖表](prediction_coverage.csv)

21 个标本的 74,888 个非空细胞全部落入冻结的每万尺度数值区间。T/NK、B/Plasma、Mono/Mac/DC 分别测得 486、485、496 个训练 HVG，分母均为 512。达到 95% 至少需要 487 个，因此 T/NK 和 B/Plasma 各两种 PCA 分类器在这两个组织均拒判。髓系 5 种方法共享 17 个标本；其余两个谱系分别由三种可预测方法共享 21 和 20 个标本。所有冻结程序的成员覆盖率为 93.1%–100%。[尺度出处](normalization_source.md) · [HVG 覆盖及中心占位](pca_measurement.csv) · [程序覆盖](../coverage/candidate_measurement.csv)

## 执行与核验

GSE245552 各谱系各抽取一个已有 OOF 袋，以开发阶段的模型和数据重算全部 15 个方法×谱系组合。与保存的 OOF logit 最大绝对差为 8.27×10⁻⁷。外部推断绑定冻结的全开发集权重、词表、训练中心和参数。按细胞 ID 对齐 scGPT 嵌入，并从同一批 ID 的程序分数重算袋均值，逐袋比对归档输入。164 上三个常规分类指标的点估计逐项与 scikit-learn 实现一致。再次按原始标本键聚合，结果与患者×组织表一致。[开发接口核验](development_interface_check.csv) · [程序袋输入](external_model_inputs.csv) · [冻结模型](../../../m5/frozen_models.csv)

外部 scGPT 在 Intel Arc A750 上按每批 4 个细胞推断，PT 和 LM 用时约 103 秒和 98 秒，两个进程峰值 RSS 均小于 3 GB；嵌入缓存共 26.1 MB。冻结参数覆盖每袋正表达条目的 94.3%–95.4%，10,693 个细胞超过 1,199 个基因的固定截断上限。[PT 编码记录](PT_embedding_manifest.csv) · [LM 编码记录](MT_embedding_manifest.csv) · [预测进程](predict_run.json)

```sh
uv run --no-project --python .venv-m4/bin/python python -B scripts/predict_gse164522.py development-check
uv run --no-project --python .venv-scgpt/bin/python python -B scripts/predict_gse164522.py encode --freeze-commit 72ee86b --tissue-code PT --device xpu --batch 4
uv run --no-project --python .venv-scgpt/bin/python python -B scripts/predict_gse164522.py encode --freeze-commit 72ee86b --tissue-code MT --device xpu --batch 4
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 uv run --no-project --python .venv-m4/bin/python python -B scripts/predict_gse164522.py predict --freeze-commit 72ee86b
```

本次结局是已取材肿瘤标本的来源。受 CD45⁺ 分选、谱系覆盖和样本规模的限制，方法之间的差异区间较宽；候选程序的患者内复现和器官背景由单独登记的检验判定。[冻结验证规则](../../../m5/validation_rules.md)
