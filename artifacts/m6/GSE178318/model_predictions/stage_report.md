# 冻结分类器（`CRLM_V1.1_245_F1`）在 GSE178318 的严格前瞻零表达暴露组织来源识别报告

> [!IMPORTANT]
> **严格前瞻零前序表达暴露（Zero Prior Expression Exposure）证明**：在 `CRLM_V1.1_245_F1` 冻结前仅审阅 `GSE178318` 公开论文与 GEO 样本/条形码元数据，未下载或触碰任何表达矩阵；全部 15 个分类器权重（`models/m4/full/` 二进制大文件位于 `.gitignore` 受控本地目录中，其逐文件 SHA256 与代码/词表/基因宇宙共 **72/72** 项哈希完整锁定于解锁前 Git 提交 `5b01467` 的 `artifacts/m5/frozen_inputs_sha256.txt` == `5ce11557b7bc565184ad7418befc08fec0fe411ab11586de7a92a8d1c7314f14`）在推理前通过 `git show 5b01467:artifacts/m5/frozen_inputs_sha256.txt` 与全量 `72/72` SHA256 逐字节校验，事后调参与重训严格为 **0**。

## 1. 核心实证结论

1. **15/15 模型全量激活，训练 HVG 外部覆盖率达 `512 / 512 = 100.0%`（零拒判、零中心替代）**：
   - 在共同核编码基因轴 `245_164_203608_shared_nuclear`（`18,733` 个核基因，与 `GSE178318` 和 `245_164_shared_nuclear` 的真实可测交集 `18,733 / 18,733` 完全一致）上，三个谱系的 `PCA_mean` 与 `PCA_attention` 均实现 **`512 / 512 (100.0%)` 精确命中**（`n_unmeasured_positions = 0`，`n_center_substituted = 0`）；19 个模型分子程序中 18 个为 `100.0%` 覆盖，仅 `B19_S4_Mac_SPP1_clean` 因 CellRanger 2.0 GRCh38 未收录克隆转录本 `AC147651.4` 为 `152 / 153 = 99.35%` 覆盖（[PCA 可测性核查](pca_measurement.csv) · [模型程序覆盖表](model_program_coverage.csv)）。
   - 对照基于真实独立分词与嵌入流水线的归档初版 `CRLM_V1_245_F1`（在原始 `33,694` 特征矩阵 `data/raw/GSE178318/` 上直接匹配原始符号）：**与 `GSE164522` 完全相同的 `4 / 15` 个模型**（`T/NK` 与 `B/Plasma` 的 `PCA_mean`、`PCA_attention`）在 `GSE178318` 上再次因原始符号覆盖率 `486 / 512 = 94.92% < 95%` 触发假性拒判；`CRLM_V1_245_F1_techfix` 恢复至 `98.83%–99.61%`（`15/15` 预测），而 `CRLM_V1.1_245_F1` 实现 **`100.0%` 全覆盖**（[三版本对比总表](three_version_gse178318_classifiers.csv)）。
2. **前瞻队列验证核心表示学结论：零样本 `scGPT` 未展现相对透明低维 `PCA` 的跨谱系增益（`T/NK` 中 `PCA` 显著优于 `scGPT_mean`）**：
   - **`Mono/Mac/DC`（6 人 / 12 部位 / 6 对 / 12 袋，入袋 `5,781` 细胞）**：**`PCA_attention` AUROC 达 `0.944`（`0.833–1.000`），`PCA_mean` 达 `0.889`（`0.778–1.000`），`scGPT_mean` 为 `0.833`（`0.722–1.000`），`scGPT_attention` 为 `0.778`（`0.556–1.000`），`program_elastic` 为 `0.750`（`0.667–1.000`）**；其中 `PCA_attention`、`PCA_mean`、`scGPT_mean` 与 `program_elastic` 四个分类器在 **`6 / 6 (100.0%)` 对同步配对患者中全部给出 `LM > PT`**（`N = 6` 下各方法配对 AUROC 差值 95% CI 跨零）。
   - **`T/NK`（6 人 / 12 部位 / 6 对 / 12 袋，入袋 `49,590` 细胞）**：**`PCA_attention` 与 `PCA_mean` AUROC 均达 `0.917`（`0.778–1.000`）**，且 **`PCA_attention` Balanced Accuracy 达 `0.833`**，两者在 **`6 / 6 (100.0%)` 对患者中全部给出 `LM > PT`**；按患者配对重抽样，**`PCA_mean` 与 `PCA_attention` 均显著优于 `scGPT_mean`（`0.583`，差值 `+0.333`，95% CI `[+0.083, +0.556]`，置信区间严格排零）**，点估计亦高于 `scGPT_attention`（`0.750`，差值 `+0.167`，95% CI `[-0.139, +0.583]`）。
   - **`B/Plasma`（6 人 / 10 部位 [`6 PT + 4 LM`] / 4 对 / 10 袋，编码 `3,967` 细胞、入袋 `3,927` 细胞；`COL12_LM=14` 与 `COL16_LM=26` 因 `<50` 细胞按规则留空）**：在全部 10 个合格部位（`method_available`）上，`PCA_mean` 与 `scGPT_mean` AUROC 均为 `0.708`，`scGPT_attention` **AUROC 为 `0.667`、Balanced Accuracy 达 `0.792`**（点估计较历史 `V1_F1` 真实子袋 AUROC `0.333`［Balanced Accuracy `0.458`］与 `Techfix` AUROC `0.333`［Balanced Accuracy `0.333`］回升；注：`N = 4` 小样本下供者配对重抽样 AUROC 差值 95% CI 分别为 `[0.000, +0.750]` 与 `[0.000, +0.722]`，仍包含 `0`），`PCA_attention` AUROC 为 `0.625` 且在 **`4 / 4 (100.0%)` 对可配对患者中全部给出 `LM > PT`**，`program_elastic` AUROC 为 `0.583`（`3 / 4` 对正向，较历史 `V1_F1`/`Techfix` 的 `0.500`［`1 / 4` 对］改善）；若限制在 4 位双侧合格配对患者的 8 个部位（`paired_donors_common_methods`），`PCA_mean`、`PCA_attention`（95% CI `0.688–1.000`）与 `scGPT_attention` AUROC 并列达 **`0.8125`**，`program_elastic` 达 **`0.750`**。
3. **术前化疗 vs 未治疗患者描述性分层结果一致（严守观测性边界，不作治疗因果推断）**：
   - 在 `treatment_naive`（`COL07, COL12, COL16`，`N = 3`，`preop_regimen = "None"`）与 `preoperative_chemotherapy`（`COL15, COL17, COL18`，`N = 3`）两个临床亚组中，`T/NK`（`PCA_mean`、`PCA_attention`）与 `Mono/Mac/DC`（`program_elastic`、`PCA_mean`、`PCA_attention`、`scGPT_mean`）均达到 **`3 / 3 (100%)` + `3 / 3 (100%)` 患者内 `LM > PT` 完全一致**；描述性上术前化疗亚组的平均 `LM - PT` logit 幅度高于未治疗亚组（如 `Mono/Mac/DC PCA_attention` `+7.557` vs `+2.185`；`T/NK PCA_attention` `+5.365` vs `+1.272`），鉴于每组仅 `N = 3` 且治疗非随机分配，仅作观察性异质性记录（[治疗分层汇总表](treatment_stratified_classifier_summary.csv)）。

## 2. 患者组织级前瞻预测总表（`scope = method_available`，10,000 次患者重抽样）

| 谱系 | 冻结分类器 | 患者/部位/配对数 | AUROC（95% bootstrap CI） | AP | Balanced accuracy | 平均 `LM − PT` logit | `LM > PT` 配对数 |
|---|---|---:|---|---:|---:|---:|---:|
| `T/NK` | `program_elastic` | 6 / 12 / 6 | 0.611（0.361–0.861） | 0.649 | 0.583 | +0.726 | 5 / 6 |
| `T/NK` | `PCA_mean` | 6 / 12 / 6 | **0.917**（0.778–1.000） | **0.915** | 0.500 | +1.813 | **6 / 6** |
| `T/NK` | `PCA_attention` | 6 / 12 / 6 | **0.917**（0.778–1.000） | **0.915** | **0.833** | +3.318 | **6 / 6** |
| `T/NK` | `scGPT_mean` | 6 / 12 / 6 | 0.583（0.444–0.778） | 0.704 | 0.417 | +0.929 | 5 / 6 |
| `T/NK` | `scGPT_attention` | 6 / 12 / 6 | 0.750（0.417–1.000） | 0.814 | 0.417 | +0.839 | 4 / 6 |
| `Mono/Mac/DC` | `program_elastic` | 6 / 12 / 6 | 0.750（0.667–1.000） | 0.751 | 0.667 | +2.423 | **6 / 6** |
| `Mono/Mac/DC` | `PCA_mean` | 6 / 12 / 6 | **0.889**（0.778–1.000） | 0.917 | **0.750** | +3.400 | **6 / 6** |
| `Mono/Mac/DC` | `PCA_attention` | 6 / 12 / 6 | **0.944**（0.833–1.000） | **0.948** | **0.750** | +4.871 | **6 / 6** |
| `Mono/Mac/DC` | `scGPT_mean` | 6 / 12 / 6 | 0.833（0.722–1.000） | 0.800 | **0.750** | +5.740 | **6 / 6** |
| `Mono/Mac/DC` | `scGPT_attention` | 6 / 12 / 6 | 0.778（0.556–1.000） | 0.766 | 0.667 | +4.543 | 5 / 6 |
| `B/Plasma` | `program_elastic` | 6 / 10 / 4 | 0.583（0.167–1.000） | 0.611 | 0.500 | +0.275 | 3 / 4 |
| `B/Plasma` | `PCA_mean` | 6 / 10 / 4 | **0.708**（0.208–1.000） | 0.715 | 0.542 | +2.412 | 3 / 4 |
| `B/Plasma` | `PCA_attention` | 6 / 10 / 4 | 0.625（0.333–0.944） | 0.542 | 0.583 | +1.711 | **4 / 4** |
| `B/Plasma` | `scGPT_mean` | 6 / 10 / 4 | **0.708**（0.333–1.000） | **0.761** | 0.500 | +1.240 | 2 / 4 |
| `B/Plasma` | `scGPT_attention` | 6 / 10 / 4 | 0.667（0.208–1.000） | 0.590 | **0.792** | +8.583 | 3 / 4 |

- [患者级完整指标](patient_metrics.csv) · [配对方法差值区间](paired_method_contrasts.csv) · [患者×组织 logits](patient_site_logits.csv) · [患者内配对差](patient_paired_logits.csv) · [治疗分层汇总](treatment_stratified_classifier_summary.csv)

## 3. 三版本（`V1_F1` / `Techfix` / `V1.1`）真实独立嵌入流水线在 GSE178318 上的平行核查

> [!NOTE]
> 三版本均采用各自真实独立的特征空间与 Intel Arc A750 XPU `scGPT` 嵌入缓存（`l_raw != l_tf != l_v11`）：
> - **`CRLM_V1_245_F1`**：读取 `data/raw/GSE178318/` 原始 `33,694` 特征计数矩阵，按 `<=256` 细胞子袋（共 `7,616` 细胞）以未映射原始无歧义符号分词编码（非零词表命中率 `94.61%–95.13%`，缓存于 `embeddings/archive/CRLM_V1_245_F1/m6/GSE178318/scgpt/`）并在全库 `log1p(CP10k)` 上做精确字符串匹配；
> - **`CRLM_V1_245_F1_techfix`**：锁定原版 M4 权重，在相同 `7,616` 个子袋细胞上经 HGNC 规范化并对齐共同核基因轴（剥离 `MT-*` / `RPS*/RPL*`，词表命中率跃升至 `99.91%–99.96%`，匹配 V1 M4 核内整数计数训练分布，缓存于 `embeddings/archive/Techfix_V1_F1/m6/GSE178318/scgpt/`）；
> - **`CRLM_V1.1_245_F1`**：前瞻冻结的洁净全袋分类头（入袋 `59,298` 细胞，`512/512 = 100.0%` HVG 覆盖，核内 `log1p(CP10k)` 预归一化 `scGPT` 嵌入）。

| 谱系 | 方法 | `CRLM_V1_245_F1`（原始 `33,694` 轴 / `7,616` 子袋细胞） | `CRLM_V1_245_F1_techfix`（HGNC+核子轴修复原权重） | `CRLM_V1.1_245_F1`（前瞻冻结全袋模型） |
|---|---|---|---|---|
| `T/NK` | `program_elastic` | 0.667（BA 0.500，6/6） | 0.667（BA 0.500，6/6） | 0.611（BA 0.583，5/6） |
| `T/NK` | `PCA_mean` | **拒判**（`486/512 = 94.92% < 95%`） | 0.889（BA 0.583，6/6，`507/512`） | **0.917**（BA 0.500，**6/6**，`512/512`） |
| `T/NK` | `PCA_attention` | **拒判**（`486/512 = 94.92% < 95%`） | **0.917**（BA 0.667，6/6，`507/512`） | **0.917**（**BA 0.833**，**6/6**，`512/512`） |
| `T/NK` | `scGPT_mean` | 0.556（BA 0.500，4/6） | 0.556（BA 0.417，4/6） | **0.583**（BA 0.417，**5/6**） |
| `T/NK` | `scGPT_attention` | 0.778（BA 0.500，5/6） | 0.722（BA 0.500，4/6） | 0.750（BA 0.417，4/6） |
| `Mono/Mac/DC` | `program_elastic` | 0.778（BA 0.750，5/6） | 0.778（BA 0.750，5/6） | 0.750（BA 0.667，**6/6**） |
| `Mono/Mac/DC` | `PCA_mean` | 0.861（BA 0.750，6/6，`496/512`） | **0.889**（BA 0.750，6/6，`510/512`） | **0.889**（BA 0.750，**6/6**，`512/512`） |
| `Mono/Mac/DC` | `PCA_attention` | 0.750（BA 0.667，4/6） | 0.750（BA 0.667，4/6） | **0.944**（BA 0.750，**6/6**，`512/512`） |
| `Mono/Mac/DC` | `scGPT_mean` | 0.833（BA 0.667，6/6） | 0.833（BA 0.750，6/6） | **0.833**（BA 0.750，**6/6**） |
| `Mono/Mac/DC` | `scGPT_attention` | 0.944（BA 0.750，5/6） | 0.917（BA 0.833，5/6） | 0.778（BA 0.667，5/6） |
| `B/Plasma` | `program_elastic` | 0.500（BA 0.500，1/4） | 0.500（BA 0.500，1/4） | **0.583**（BA 0.500，**3/4**；配对集 `0.750`） |
| `B/Plasma` | `PCA_mean` | **拒判**（`486/512 = 94.92% < 95%`） | **0.875**（BA 0.542，4/4，`506/512`） | 0.708（BA 0.542，3/4；配对集 `0.8125`） |
| `B/Plasma` | `PCA_attention` | **拒判**（`486/512 = 94.92% < 95%`） | 0.792（BA 0.542，4/4，`506/512`） | 0.625（BA 0.583，**4/4**；配对集 `0.8125`） |
| `B/Plasma` | `scGPT_mean` | 0.750（BA 0.500，2/4） | 0.708（BA 0.500，2/4） | 0.708（BA 0.500，2/4） |
| `B/Plasma` | `scGPT_attention` | 0.333（BA 0.458，2/4） | 0.333（BA 0.333，2/4） | **0.667**（**BA 0.792**，**3/4**；配对集 `0.8125`） |

## 4. 标本边界、编码细胞 vs 入袋细胞区分与 Intel Arc A750 XPU 记录

1. **PBMC 物理阻断与编码细胞（`59,338`）vs 入袋细胞（`59,298`）精确区分**：
   - 仅纳入 6 位同步 IV 期 MSS 结直肠癌肝转移患者的 `12` 个实体瘤文库（`6 PT + 6 LM`：`GSM5387660–663`、`GSM5387665–670`、`GSM5387672–673`），3 个外周血文库（`GSM5387664`、`GSM5387671`、`GSM5387674`）在入口清单与目录扫描双重断言下 `100%` 阻断；
   - **全量编码细胞数 = `59,338`**：12 个实体瘤文库中全部通过 M2 QC 的核心三谱系细胞（`T/NK = 49,590`、`Mono/Mac/DC = 5,781`、`B/Plasma = 3,967`）均已编码缓存于 `embeddings/m6/GSE178318/scgpt/GSM5387660–673.npz`；
   - **合格全袋数 = `34 / 36`、入袋细胞数 = `59,298`**：按 `>= 50` 个 QC 保留细胞门槛建袋，`T/NK`（`12/12` 袋，`49,590` 细胞）与 `Mono/Mac/DC`（`12/12` 袋，`5,781` 细胞）全数入袋，`B/Plasma` 共 `10/12` 袋（`3,927` 细胞）入袋，仅 `GSM5387663`（`COL12_LM`，`14` 细胞）与 `GSM5387668`（`COL16_LM`，`26` 细胞）共 `40` 个细胞因 `< 50` 细胞未入袋（`59,338 - 40 = 59,298`）；`180` 次模型×标本判定中 **`170` 次输出预测，`10` 次仅因 `< 50` 细胞留空，`0` 次因特征缺失拒判**（[预测覆盖表](prediction_coverage.csv) · [袋级预测](bag_predictions.csv) · [细胞名单](selected_cells.csv.gz) · [外部程序袋均值](external_program_bag_scores.csv)）。
2. **Intel Arc A750 XPU 编码与资源实测**：
   - `V1.1` 全袋核轴编码：在 Intel Arc A750 XPU（`batch=16`，`OMP/OPENBLAS/MKL_NUM_THREADS=4`）上编码 `59,338` 个细胞，模型纯前向总耗时 `726.39s`（吞吐率 **`81.7 cells/s`**，总文件体积 `113.9 MB`），非零表达词表命中率 **`99.90%–99.96%`**（[编码清单](embedding_manifest.csv)）；
   - 历史三版本对照子袋编码：在相同 A750 XPU 上独立编码 `7,616` 个子袋细胞的原版 `V1_F1` 原始 `33,694` 轴嵌入（`94.97s`）与 `Techfix_V1_F1` 核轴整数计数嵌入（`92.30s`），连同缓存校验与流式提取总耗时 `225.72s`（峰值 RSS `3.01 GiB`，[encode_run.json](encode_run.json)）；
   - 全袋推理、10,000 次患者级重抽样与真实三版本对比总耗时 `29.26s`（[predict_run.json](predict_run.json) · [predict.log](predict.log)）。

```bash
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
.venv-scgpt/bin/python -B scripts/predict_gse178318.py --freeze-commit 5b01467
sha256sum -c --quiet artifacts/m6/GSE178318/model_predictions/output_sha256.txt
```
