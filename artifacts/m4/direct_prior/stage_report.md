# M4 · 层级 MIL（HA / CTA）与变分混合效应 MixMIL 直接对照报告（Track 2 · `CRLM_V1.1`）

> [!NOTE]
> 初版 `CRLM_V1_245_F1`（单谱系、256 细胞截断、未调参固定 100 步 HA/CTA 与含 `[::2]` 切片缺陷的 MixMIL）已完整封存于 `artifacts/archive/CRLM_V1_245_F1/artifacts/m4/direct_prior/` 与 `models/archive/CRLM_V1_245_F1/m4_direct/`。本报告记录 `V1.1` 方法学重构版在全量合格细胞袋、洁净分层细类及跨三谱系联合全标本袋（`All3Lineages`）上的完整对照结果。

## 1. 核心结论

在严格一致的外层五折留患者划分、共同核编码基因轴（`245_164_shared_nuclear`，`18,815` 个基因）与全量合格细胞袋（`70,628` 个细胞）上，即使补齐了 Do & Lähdesmäki (2025) 原设定的**跨三谱系全标本 11 细类层级结构（`All3Lineages`）**、折内 3 折患者互斥调参（`dropout`、`weight_decay`、`steps`）与修复后的全协方差变分 `MixMIL`（移除 `[::2]` 奇偶丢袋并启用 `1/sqrt(Q)` 注意力缩放）：
- **贝叶斯收缩 `MixMIL` 与透明均值池化表现稳健领先**：在 `All3Lineages`（16 位患者、29 个患者部位、13 对 PT–LM）上，`PCA32_MixMIL` 折外患者组织级 AUROC 达 **`1.000`**（相对同表示 `PCA_mean = 0.986` 差值 `+0.014`，95% CI `[0.000, +0.048]`），髓系单谱系 `PCA32_MixMIL = 0.994`（相对 `PCA_mean = 0.988` 差值 `+0.006`，95% CI `[0.000, +0.028]`）；`scGPT32_MixMIL` 在 `All3Lineages` 达 **`0.957`**（13/13 配对患者全员 `LM > PT`），在髓系达 **`0.941`**。
- **内层调参与洁净细类显著改善 `HA` / `CTA`，但在 `N=15–16` 小患者样本量下仍不优于均值池化**：例如 `T/NK` 的 `PCA32_CTA` 由初版 `0.759` 升至 **`0.938`**、`scGPT512_CTA` 由 `0.605` 升至 **`0.805`**，`B/Plasma` 的 `PCA32_CTA` 由 `0.658` 升至 **`0.767`**；但在全部四项任务与三种表示下，两层非线性层级注意力网络（`1,122–16,515` 参数）均未显著超越对应线性均值头（`33–513` 参数），符合 `AGENTS.md` 第 7 节保留简约基线的原则。

## 2. 四项任务折外患者组织级 AUROC 总表（主种子 `20261121`）

| 输入表示与模型 | `T/NK` (16人/28部位/12对) | `Mono/Mac/DC` (15人/26部位/11对) | `B/Plasma` (15人/23部位/8对) | `All3Lineages` 联合标本 (16人/29部位/13对) |
|---|---:|---:|---:|---:|
| `PCA_mean`（匹配均值头） | **0.995** | 0.988 | **0.808** | 0.986 |
| `PCA_attention`（透明线性注意力） | 0.903 | 0.982 | 0.658 | 0.986 |
| `PCA32_CTA`（折内调参细类注意力） | 0.938 | 0.882 | 0.767 | 0.833 |
| `PCA32_HA`（折内调参双层注意力） | 0.944 | 0.888 | 0.600 | 0.900 |
| `PCA32_MixMIL`（修复版高斯变分 GLMM） | 0.985 | **0.994** | 0.783 | **1.000** |
| `scGPT32_mean`（匹配均值头） | **0.908** | **0.953** | **0.833** | **0.971** |
| `scGPT32_attention`（透明线性注意力） | 0.523 | 0.811 | 0.450 | 0.881 |
| `scGPT32_CTA`（折内调参细类注意力） | 0.774 | 0.746 | 0.600 | 0.738 |
| `scGPT32_HA`（折内调参双层注意力） | 0.738 | 0.757 | 0.600 | 0.724 |
| `scGPT32_MixMIL`（修复版高斯变分 GLMM） | 0.887 | 0.941 | 0.725 | 0.957 |
| `scGPT_mean`（原生 512D 均值头） | **0.892** | 0.893 | **0.817** | **0.981** |
| `scGPT_attention`（原生 512D 线性注意力） | 0.644* | **0.911** | 0.300 | 0.910 |
| `scGPT512_CTA`（折内调参细类注意力） | 0.805 | 0.811 | 0.575 | 0.771 |
| `scGPT512_HA`（折内调参双层注意力） | 0.759 | 0.882 | 0.658 | 0.790 |

*\*注：`T/NK` 的 512 维无投影 `scGPT_attention`（1,025 参数）在第 0 折达 500 步 L-BFGS-B 迭代上限而按预定规则保留未收敛状态（其余 4 折 13 人/23 部位有效）；所有 `PCA32`、`scGPT32` 及 `HA/CTA/MixMIL` 分支均为 100% 完整患者覆盖。*

- [患者级完整指标](patient_metrics.csv) · [相对同表示均值头的配对区间](matched_mean_comparisons.csv) · [患者组织 logits](patient_site_logits.csv)

## 3. 方法学保真度、防泄漏内层调参与上游数值等价核验

1. **仓内源码锁定（零 `~/.cache` 隐式依赖）**：直接绑定并经 [PROVENANCE.json](../../../references/code/PROVENANCE.json) 校验 SHA256 的 `references/code/hier-mil`（commit `2fb155840394961595755bede506c4b54c0942ae`）与 `references/code/mixmil`（commit `16d92349ccf19b73c1dccd167934fcc3202b6be4`）。
2. **双层生物层级与严格折内 3 折调参**：
   - 单谱系采用 M2 零 `unresolved` 细类标签；`All3Lineages` 合并全部 11 类 `lineage::fine_subtype` 免疫亚群；
   - `HA` 与 `CTA` 严格在每折 `train_donors` 内部按配对状态构造 3 折患者互斥验证集（每个内层训练子集独立重拟合 HVG/PCA 或 `scGPT32` 投影与标准化），在 8 组配置 `dropout in {0.0, 0.3} x weight_decay in {1e-3, 1e-2} x steps in {100, 250}` 上选优后于外层训练折重拟合（调参轨迹完整记录于 `PCA32/inner_cv_tuning.csv`、`scGPT32/inner_cv_tuning.csv`、`scGPT512/inner_cv_tuning.csv`）。
3. **MixMIL 修复与上游张量核等价性**：
   - `src/crlm/prior_mil.py` 已修复 `mean_glmm_init` 误切片 `[::2]` 并开启 `scale_by_sqrt_q=True`；
   - `check_upstream.py` 实测 `CTA`、`HA` 与上游 `hier-mil/model.py` 最大绝对误差均为 **`0.0`**，修复版 `MixedMIL` 与上游 `MixMIL._calc_bag_emb_effect_tensor`（含 `1/sqrt(Q)` 缩放）最大绝对误差为 **`2.38 × 10⁻⁷`**（[upstream_equivalence.json](upstream_equivalence.json)）；
   - 全部折外测试袋逐袋通过单标本独立推理一致性断言、注意力权重和 `== 1` 断言及带符号细胞贡献加截距精确重构袋 logit 断言。
4. **运行资源**：`PCA32`、`scGPT32`、`scGPT512` 分别执行 `85`、`120`、`60` 次外层模型拟合（另含 `960` 组内层 CV 评估，在 Intel Arc A750 XPU 上加速 `HierarchicalMIL`），耗时分别为 `895.8s`、`848.2s`、`1170.8s`，峰值 RSS 分别为 `8.01 GiB`、`4.98 GiB`、`5.33 GiB`。

## 4. 可复现执行命令

```bash
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
for rep in PCA32 scGPT32 scGPT512; do
    .venv-scgpt/bin/python -Bu scripts/benchmark_prior_mil.py --representation "$rep" --device xpu
done
.venv-scgpt/bin/python -Bu artifacts/m4/direct_prior/check_upstream.py
.venv-m4/bin/python -Bu scripts/summarize_benchmarks.py
```
