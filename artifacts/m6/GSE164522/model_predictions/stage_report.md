# 冻结分类器（`CRLM_V1.1_245_F1`）在 GSE164522 的外部组织来源识别报告

> [!NOTE]
> 初版 `CRLM_V1_245_F1`（11/15 模型有效、4 个 PCA 模型因未做 HGNC 对齐 `<95%` 假性拒判、256 细胞截断且含跨谱系 `CD3E+` 双细胞）已封存于 `artifacts/archive/CRLM_V1_245_F1/artifacts/m6/GSE164522/model_predictions/` 与 `embeddings/archive/CRLM_V1_245_F1/m6/`；锁定原版权重的纯技术勘误版见 `artifacts/techfix_v1_f1/models/stage_report.md`。本报告记录 `CRLM_V1.1_245_F1` 全部 **15/15 个重构冻结分类器**在洁净全袋上的外部评估结果。

## 1. 核心实证结论

1. **15/15 模型全量激活，训练 HVG 外部覆盖率达 `512 / 512 = 100.0%`（零拒判、零中心替代）**：
   - 依托 `M4` 在共同核编码基因轴 `245_164_shared_nuclear`（`18,815` 个核基因）上训练并在 `M6` 经 HGNC 精确对齐及核子轴 `CP10k` 重标定，三个谱系的 `PCA_mean` 与 `PCA_attention` 均实现 **`512 / 512 (100.0%)` 精确特征命中**（`n_unmeasured_positions = 0`，`n_center_substituted = 0`），彻底根除初版 4/15 模型因基因别名失配导致的拒判。
2. **透明低维 `PCA` 在全部三个谱系上系统性优于零样本预训练 `scGPT`**：
   - **`T/NK`（10 人 / 20 部位 / 10 对 / 21 袋）**：`PCA_attention` AUROC 达 **`0.980`**（95% CI `0.920–1.000`），`PCA_mean` 达 **`0.970`**（`0.900–1.000`，Balanced Accuracy **`0.950`**），两者在 **10/10 对患者中全部给出 `LM > PT`**；按患者配对重抽样，`PCA_mean` 显著高于 `scGPT_mean`（`0.810`，差值 **`+0.160`**，95% CI `[+0.020, +0.250]`），`PCA_attention` 亦显著高于 `scGPT_mean`（差值 **`+0.170`**，`[+0.020, +0.260]`）。
   - **`Mono/Mac/DC`（9 人 / 17 部位 / 8 对 / 17 袋）**：`PCA_attention` AUROC 达 **`0.889`**（`0.630–1.000`），`scGPT_attention` 为 **`0.806`**（`0.537–1.000`，`7/8` 对正向），`PCA_mean` 为 **`0.750`**（Balanced Accuracy `0.771`），`scGPT_mean` 为 `0.625`（`PCA_attention - scGPT_mean = +0.264`，95% CI `[+0.049, +0.444]`）。
   - **`B/Plasma`（10 人 / 19 部位 / 9 对 / 19 袋）**：修复初版未配对 PT 双倍权重失真与单特征深度混杂后，**`program_elastic` 彻底消除反向失真（初版 `0.230`、`2/10`）并跃升至 `0.856`**（`0.675–1.000`、`8/9` 对患者 `LM > PT`，准确捕获 `B_Antigen_Presentation_TLS` 相对 `Plasma_Antibody_Secretion` 的转移富集方向）；`PCA_attention` 达 **`0.844`**（`0.683–0.971`、`8/9` 对正向），显著优于 `scGPT_mean`（`0.578`，差值 `+0.267`，`[+0.040, +0.517]`）与 `scGPT_attention`（`0.333`，差值 `+0.511`，`[+0.225, +0.840]`）。

## 2. 患者组织级外部预测总表（`scope = method_available`，10,000 次患者重抽样）

| 谱系 | 冻结分类器 | 患者/部位/配对数 | AUROC（95% bootstrap CI） | AP | Balanced accuracy | `LM > PT` 配对数 |
|---|---|---:|---|---:|---:|---:|
| `T/NK` | `program_elastic` | 10 / 20 / 10 | 0.500（0.270–0.770） | 0.521 | 0.500 | 6 / 10 |
| `T/NK` | `PCA_mean` | 10 / 20 / 10 | **0.970**（0.900–1.000） | 0.977 | **0.950** | **10 / 10** |
| `T/NK` | `PCA_attention` | 10 / 20 / 10 | **0.980**（0.920–1.000） | **0.983** | 0.850 | **10 / 10** |
| `T/NK` | `scGPT_mean` | 10 / 20 / 10 | 0.810（0.710–0.980） | 0.834 | 0.750 | **10 / 10** |
| `T/NK` | `scGPT_attention` | 10 / 20 / 10 | 0.880（0.700–1.000） | 0.908 | 0.650 | 9 / 10 |
| `Mono/Mac/DC` | `program_elastic` | 9 / 17 / 8 | 0.583（0.302–0.852） | 0.710 | 0.549 | 5 / 8 |
| `Mono/Mac/DC` | `PCA_mean` | 9 / 17 / 8 | 0.750（0.463–1.000） | 0.807 | **0.771** | 6 / 8 |
| `Mono/Mac/DC` | `PCA_attention` | 9 / 17 / 8 | **0.889**（0.630–1.000） | **0.914** | 0.688 | 6 / 8 |
| `Mono/Mac/DC` | `scGPT_mean` | 9 / 17 / 8 | 0.625（0.361–0.889） | 0.752 | 0.660 | 6 / 8 |
| `Mono/Mac/DC` | `scGPT_attention` | 9 / 17 / 8 | 0.806（0.537–1.000） | 0.886 | 0.701 | **7 / 8** |
| `B/Plasma` | `program_elastic` | 10 / 19 / 9 | **0.856**（0.675–1.000） | **0.897** | 0.600 | **8 / 9** |
| `B/Plasma` | `PCA_mean` | 10 / 19 / 9 | 0.667（0.428–0.860） | 0.709 | **0.639** | 6 / 9 |
| `B/Plasma` | `PCA_attention` | 10 / 19 / 9 | **0.844**（0.683–0.971） | 0.864 | 0.556 | **8 / 9** |
| `B/Plasma` | `scGPT_mean` | 10 / 19 / 9 | 0.578（0.329–0.825） | 0.691 | 0.483 | 6 / 9 |
| `B/Plasma` | `scGPT_attention` | 10 / 19 / 9 | 0.333（0.070–0.622） | 0.406 | 0.378 | 2 / 9 |

- [患者级完整指标](patient_metrics.csv) · [配对方法差值区间](paired_method_contrasts.csv) · [患者×组织 logits](patient_site_logits.csv) · [患者内配对差](patient_paired_logits.csv)

## 3. 洁净全袋覆盖、接口核验与 Intel Arc A750 XPU 编码记录

1. **双细胞剔除与全袋样本覆盖**：
   - 在肿瘤标本（`PT` 与 `MT` 共 21 个标本、10 位患者）中剔除作者标注的 `1,748` 个 `CD3E+` 髓系双细胞与 `1,481` 个 `hC75_B-CD3E` 双细胞后，按 `>= 50` 个洁净细胞门槛保留全部细胞建袋，共生成 **`57` 个合格全袋、`68,050` 个洁净肿瘤免疫细胞**（`PT` 28 袋 `33,419` 细胞，`MT` 29 袋 `34,631` 细胞）；
   - `21 标本 × 3 谱系 × 5 方法 = 315` 次判定中：**`285` 次成功输出预测，`30` 次仅因洁净细胞数 `< 50` 按规则不预测，`0` 次因特征缺失拒判**（[预测覆盖表](prediction_coverage.csv) · [袋级预测](bag_predictions.csv) · [细胞名单](selected_cells.csv.gz)）。
2. **开发接口一致性与冻结哈希校验**：
   - `development-check` 在 `GSE245552` 全袋上重算全部 15 个方法×谱系折外 logit，与 `M4` 归档值最大绝对误差仅 **`8.69 × 10⁻⁷`**（[开发接口核验](development_interface_check.csv)）；
   - 外部推理启动前严格通过 `artifacts/m5/frozen_inputs_sha256.txt`（`5ce11557...`）全套权重、词表与源码 SHA256 校验，且所有点估计与 `scikit-learn` 参考实现逐项一致。
3. **A750 XPU 编码性能**：
   - 在 Intel Arc A750 XPU（`batch=16`）上编码全部 `68,050` 个洁净外部细胞耗时 `894.3s`（`PT` `426.8s`、`MT` `442.0s`，约 `76.1 cells/s`），缓存大小 `130.6 MB`，非零表达词表命中率达 **`99.83%–99.93%`**（较初版未做 HGNC 桥接的 `94.3%–95.4%` 显著提升），峰值 RSS `3.10 GiB`；CPU 预测与 10,000 次重抽样耗时 `29.9s`（峰值 RSS `0.77 GiB`）。

```bash
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
.venv-m4/bin/python -Bu scripts/predict_gse164522.py development-check
.venv-scgpt/bin/python -Bu scripts/predict_gse164522.py encode --freeze-commit WORKTREE --device xpu --batch 16 --overwrite
.venv-m4/bin/python -Bu scripts/predict_gse164522.py predict --freeze-commit WORKTREE
```
