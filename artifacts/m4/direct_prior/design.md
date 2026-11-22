# 直接前作（HA / CTA / MixMIL）对照设计（Track 2 · `CRLM_V1.1`）

## 1. 仓内锁定源码与四类评估任务

- **零仓外隐式依赖**：全部上游参考实现封存于仓内 `references/code/hier-mil`（commit `2fb155840394961595755bede506c4b54c0942ae`，Do & Lähdesmäki 2025）与 `references/code/mixmil`（commit `16d92349ccf19b73c1dccd167934fcc3202b6be4`，Engelmann et al. 2024），运行时对照 `references/code/PROVENANCE.json` 校验关键源码 SHA256。
- **单谱系 + 跨谱系联合样本级（`All3Lineages`）双层设计**：
  1. **单谱系全袋任务**（`T/NK`、`Mono/Mac/DC`、`B/Plasma`）：使用 M2 洁净细类标签（`CD8 T cell`、`CD4 T cell`、`Double_neg_or_unsplit_T`、`NK cell`；`Macrophage`、`Monocyte`、`cDC`、`mature DC`、`pDC`；`B cell`、`plasma cell`）作为袋内子层级；
  2. **跨三谱系全标本联合任务（`All3Lineages`）**：合并同一肿瘤文库的全部 `T/NK + Mono/Mac/DC + B/Plasma` 合格细胞（33 个文库、70,628 个细胞、11 个 `lineage::fine_subtype` 免疫细类），严格复现 Do & Lähdesmäki (2025)“细胞 → 细胞类型 → 标本表型”的多谱系层级 MIL 设定，并同步运行同输入的 `PCA_mean/attention`、`scGPT32_mean/attention`、`scGPT_mean/attention` 配对基线。

## 2. 折内嵌套内层 CV 调参与 MixMIL 缺陷修复

- **HA / CTA 折内 3 折患者分组调参（零外层泄漏）**：
  - 严格在每个外层训练折 `train_donors` 内部构造按配对状态分层的 3 折患者互斥内层划分，并在每个内层训练子集上独立重拟合特征投影与标准化（严禁外层测试折或内层验证折参与缩放）；
  - 在 `dropout in {0.0, 0.3} x weight_decay in {1e-3, 1e-2} x steps in {100, 250}`（共 8 组配置，完整记录于 `inner_cv_tuning.csv`）上按内层折外患者组织级 AUROC（平局偏向更强正则化）选定超参数后，在外层 `train_donors` 重拟合。
- **MixMIL 修复版（`MixedMIL`）**：
  - 修复初版 `mean_glmm_init` 在 `n_trials=1` 下误照搬上游 `n_trials=2` 的 `[::2]` 切片（恢复全量训练袋方差比估计），并补齐原论文 `1 / sqrt(Q)` 注意力特征缩放（`scale_by_sqrt_q=True`）；
  - 在 `PCA32` 与 `scGPT32`（`64 x 64` 全协方差高斯变分后验，500 步 MC ELBO + 患者配对项）上运行单谱系与 `All3Lineages` 四项任务，并由 `check_upstream.py` 核验与上游张量核的逐元素等价性（`< 2e-6`）。
