# M6 · 独立队列与四组织背景

两个候选在主要外部队列 GSE164522 的平均值和患者多数方向均与发现集一致。两项均未达到冻结的主要检验标准。UPenn 次要来源的独立患者差值分向两侧。验证版本为 `CRLM_V1_245_F1`，完整八项登记中六项有实测结果，另外两项随 GSE178318 使用条件保留缺失。[登记检验](registered_tests/registered_test_results.csv) · [检验家族](registered_tests/registered_families.csv)

| GSE164522 患者肿瘤差 | 独立患者 | 同向 | 平均 LM−PT | 患者 bootstrap 95% 区间 | 双侧 sign p | 族内 Holm p |
|---|---:|---:|---:|---:|---|---:|
| Mono/Mac/DC Mac_SPP1 | 8 | 6 上升 | +0.00981 | −0.00024 至 +0.02100 | 0.2891 | 0.2891 |
| B/Plasma UPR | 10 | 8 下降 | −0.01577 | −0.02544 至 −0.00541 | 0.1094 | 0.2188 |

符号检验评估患者方向，bootstrap 区间描述均值的不确定性。两者使用不同的统计量。原始表达来自四张连续值矩阵，核对 163,347 个细胞与原 metadata，42 个标本、10 位患者。共同排序轴 24,662 个基因，两个程序的覆盖分别为 149/160 和 106/107。主要检验用每标本每谱系至少 50 个细胞。[基因覆盖](GSE164522/coverage/candidate_measurement.csv) · [患者表](patient_contrasts.csv)

| 同一批四组织患者 | 四组织 N | LM−PT | NL−NC | 交互 I | I 同向 | 背景族 Holm p |
|---|---:|---:|---:|---:|---:|---:|
| Mac_SPP1 | 5 | +0.01637 | +0.00104 | +0.01533 | 5/5 | 0.125 |
| UPR | 10 | −0.01577 | −0.01308 | −0.00269 | 6/10 | 0.754 |

UPR 的正常器官差与肿瘤差幅度相近。把 B 和 Plasma 分开，两个亚类均有支持时的患者平均组成项为 −0.0111、类型内状态项为 −0.0047；在每档至少一个细胞的深度分解中，UPR 的档内项为 −0.0146。髓系在排除作者 CD3E 命名亚组后，8 人中有 5 人同向，均值 +0.0076。组成、深度与外部亚型界定均影响终点解释。[敏感性](GSE164522/patient_effects/sensitivity_summary.csv) · [B 与浆细胞分解](GSE164522/patient_effects/patient_B_plasma_decomposition.csv) · [检出深度分解](GSE164522/patient_effects/patient_detection_depth_decomposition.csv)

GSE203608 髓系 4 人为 2 升 2 降，均值 −0.00499；UPR 的两位合格患者一升一降。患者 49 的临床时序敏感性保留在次要报告中。GSE178318 作者在 GEO 要求与 PI 确认使用条件，本版本该来源未测，两项原始 p 为 NA，仅为保留家族规模在 Holm 输入中记 1。[UPenn 结果](GSE203608/stage_report.md) · [GSE178318 使用条件](../../starter/audit_data/geo/GSE178318_metadata_only.soft)

冻结分类器在 164 的 21 个肿瘤标本上进行迁移。15 个头中有 11 个满足测量条件并产生 208 条袋级预测；T/NK 和 B/Plasma 的训练 PCA 基因覆盖分别为 486/512 和 485/512，低于冻结的 95% 门槛。髓系为 496/512。每个谱系各方法在相同可测患者上的比较另表保存。[覆盖及拒判](GSE164522/model_predictions/prediction_coverage.csv) · [患者指标与区间](GSE164522/model_predictions/patient_metrics.csv)

使用 `scripts/assemble_validation.py` 汇总各队列真实患者表，逐患者核对四组织交互恒等式，再执行 starter 冻结家族校正。八项登记全部保留，十四项未登记的描述性对比另表归档。[患者效应](effects_summary.csv) · [校正审计](registered_tests/registered_test_audit.json) · [描述性对比](registered_tests/unplanned_exploratory_results.csv)
