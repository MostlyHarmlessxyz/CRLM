# 发现集候选草案

本表由 M3 登记的四个完整文献程序产生。患者内差为 LM 减 PT，单位为 0–1 秩分数。

| 谱系 | 程序 | 完整患者对 | 平均差 | 同向患者 | 留一范围 | 短名单 |
|---|---|---:|---:|---:|---|---|
| Mono/Mac/DC | B19_S4_Mac_SPP1 | 11 | 0.01974 | 11/11 | 0.01530 至 0.02150 | 入选 |
| B/Plasma | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 8 | -0.01219 | 7/8 | -0.01909 至 -0.00998 | 入选 |
| Mono/Mac/DC | B19_S4_Mac_CXCL10 | 11 | -0.00571 | 7/11 | -0.01021 至 -0.00239 | 同向患者少于三分之二 |
| T/NK | HALLMARK_INTERFERON_GAMMA_RESPONSE | 12 | 0.00138 | 6/12 | -0.00094 至 0.00384 | 同向患者少于三分之二、留一患者后方向变化 |

每个候选的门槛 20/100、宽松 QC、去共享成员、去肝/血红蛋白及 DRB3/4/5、检测率匹配随机集合的结果见[敏感性结果](patient_effects/sensitivity_summary.csv)。同一父谱系内的[患者分解](composition_results/patient_subtype_decomposition.csv)列出共同支持质量。冻前方向按患者均差确定，模型增量归 M4，确认家族归 M5。

测序深度、细胞组成和伪批量的独立分析见[阶段报告](stage_report.md)。
