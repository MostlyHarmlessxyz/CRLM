# GSE203608 · UPenn 独立体内次要来源（V1.1 `CRLM_V1.1_245_F1`）

在 M5 原子冻结（`CRLM_V1.1_245_F1`）后，基于官方 HGNC 规范化与 Ensembl `gene_id` 消歧（`33,694` 行源特征消歧得到 `33,658` 个唯一可测规范基因，投影至三队列共同核基因轴 `245_164_203608_shared_nuclear` 共 `18,733` 个核基因）、严格细胞级质量控制（`Cell_type != "NA"` 且 `MitoFraction <= 0.20`，从 `59,661` 个体内细胞中剔除 `12,136` 个死亡/未标注高线粒体细胞，保留 **47,525** 个洁净体内细胞），使用开发集冻结 `matched_controls.csv`（零外部重拟合）完成次要方向与精度评估。

## 预登记次要检验与焦点细类结果

- **`Mono/Mac/DC::Core_Lipid_SAM_R`（预期 `positive`，实测成员 `16/16 = 100.0%`）**：
  - 宽谱系 `Mono/Mac/DC`（`N=4` 对解剖 PT–LM）：均值 **`+0.01713`**（2 正 / 2 负，双侧精确 `p = 1.000`，Holm `p = 1.000`）；剔除原作者图示临床 Stage III 而标本表含 `LiverMet` 的 `Patient 49` 后（`N=3`），均值升至 **`+0.07141`（2 正 / 1 负）**。
  - **焦点细类纯 `Macrophage`（剥离 DC 稀释效应后，`N=4`）**：在洁净作者标注下均值升至 **`+0.03051`（3 正 / 1 负，75% 同向）**，在跨队列统一 `annotation_m1.classify` 规则下同样为 **`3 正 / 1 负`（`+0.01134`；`monocyte/macrophage` `+0.01593`，3 正 / 1 负）**！
- **`T/NK::B07_Th17_TWEAK`（预期 `negative`，实测成员 `14/14 = 100.0%`）**：
  - 宽谱系 `T/NK`：主门槛 `>=50` 细胞下仅 `N=2` 对患者满足双侧支持，均值 **`−0.00850`**（1 正 / 1 负，`p = 1.000`，Holm `p = 1.000`）；放宽至 `>=20` 细胞门槛（`N=3`）时均值为 **`−0.01433`（1 正 / 2 负）**，在跨队列统一 `annotation_m1` 规则 `>=20` 门槛下（`N=3`）达 **`0 正 / 3 负`（均值 `−0.02169`，`raw p = 0.250`）**，并在原本作者未拆分 CD4/CD8 的 `GSE203608` 中成功解析出纯 `CD4 T cell` 对比（`>=50`: `N=1, −0.02738`; `>=20`: `N=2, −0.00869`）。

[登记次要检验](registered_tests/registered_results.csv) · [逐患者数据](patient_effects/patient_contrasts.csv) · [剔除 Patient 49](patient_effects/exclude_patient49.csv) · [敏感性全表](patient_effects/sensitivity_summary.csv) · [焦点细类表](patient_effects/focal_subtype_effects_summary.csv)

## 数据核验与执行资源

输入为原版 `126,014` 个细胞的整数计数矩阵（`315,102,049` 个非零坐标全部校验通过）及完整原版细胞注释（与 7 列历史身份表逐行全等且完整保留 `MitoFraction`）。24 个体内标本覆盖 15 位患者代码（4 对解剖 PT–LM）。双免疫组成分母与作者/规则注释对照表一并归档于 `coverage/`。全流程耗时 63.51 秒，峰值 RSS 4.64 GiB。[细胞 QC 审计](coverage/cell_qc_audit.csv) · [基因覆盖](coverage/measurement_coverage.csv) · [双分母组成](coverage/sample_immune_composition.csv) · [资源记录](resource_profile.json)
