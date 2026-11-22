# M4 开发集模型比较约定（Track 2 · `CRLM_V1.1_AUDIT_REMEDIATED`）

> [!NOTE]
> 初版 `CRLM_V1_245_F1` 的 M4 产物、模型与 256 细胞抽样嵌入已完整封存于 `artifacts/archive/CRLM_V1_245_F1/artifacts/m4/`、`models/archive/CRLM_V1_245_F1/` 与 `embeddings/archive/CRLM_V1_245_F1/`。本文档记录 `V1.1` 方法学重构版的开发集建模规约。

## 1. 样本袋、共同核基因轴与零正常组织掺入

- **独立单位与外层分折**：沿用 M0 预登记的 5 个患者折（`data/manifests/m0/donor_folds_pre_qc.json`），独立单位严格为核实后的患者 `donor_uid`。同患者所有病灶与文库绑定同一训练/测试折。
- **肿瘤分类边界**：仅 `GSE245552` 的 `PT` 与 `LM` 肿瘤文库可进入分类器训练与评估（`src/crlm/model_bags.py` 对 `NC`/`NL` 正常组织强制抛出 `ValueError`）。
- **全袋保留与抽样方差审计**：每个肿瘤文库的 `T/NK`、`Mono/Mac/DC`、`B/Plasma` 在 `qc_retained` 细胞数 `>= 50` 时保留**全部合格细胞**（`max_cells=None`，共 `88` 个单谱系袋、`70,388` 个细胞；另建 `All3Lineages` 三谱系联合标本袋共 `33` 个袋、`70,628` 个细胞），彻底消除初版单次 256 细胞截断损失；同时附 `256 x 10` 次独立确定性子袋抽样敏感性（`subbag_256x10_sensitivity.csv`）量化抽样方差。
- **跨平台共同核编码基因轴**：所有特征统一对齐至 HGNC 规范化的 `245_164_shared_nuclear`（`18,815` 个核编码基因，完全剔除正常组织特有的 `353` 个线粒体/核糖体偏倚基因），从构造上保证训练出的全部 512 个 PCA HVG 在外部 `GSE164522` 上达到 `512/512 = 100.0%` 覆盖（消除初版 4/15 模型假性拒判）。

## 2. 折内严格防泄漏特征工程与线性基线

- **程序与组成 Elastic Net 基线（消除未配对 PT 双倍权重失真）**：
  - `program_elastic`：采用 M3 预定义的 19 个完整程序（髓系 10、T/NK 6、B/Plasma 3）；**严格在每个外层训练折 `train_donors` 内部**计算基因检出率分桶（20 桶）与 `n_ctrl=10` 同桶非成员背景对照，再对训练与测试袋统一计算深度去偏 `debiased_ucell_scores`（全集迁移头在全开发集锁定 `program_defs.json` 与 `matched_controls.csv` 供 M6 冻结调用）。
  - `composition_elastic`：采用标本全量 `qc_retained` 细类比例（`specimen_immune_composition`，零抽样噪声）。
  - `QC_elastic`：仅用各袋 `log1p(n_counts_visible)` 与 `log1p(n_genes_visible)` 的四分位数检验测序深度混杂基线。
  - **类别均衡患者部位权重 + 配对排序项**：令训练折内 `PT` 患者集与 `LM` 患者集总权重严格各占 `0.5`（同患者同组织内多病灶等权，消除初版未配对 PT 患者获得 2 倍权重导致的 `B/Plasma` 系数符号反转），并联合 `lambda_pair=0.25` 患者内 `LM - PT` 配对 softplus 损失，以盒约束 `L-BFGS-B` 精确求解 `L1 + L2`（`alpha=0.05, l1_ratio=0.5`）。
- **PCA32 与折内实质污染掩码（`PCA_remove_organ_markers`）**：
  - 在 `18,815` 核基因轴上逐细胞计算 `log1p(CP10k)`，严格在每折 `train_donors` 内选方差最高的 512 个基因并拟合 32 维中心化 PCA。
  - `PCA_remove_organ_markers`：严格仅读取 `train_donors` 内的 `NL` 与 `NC` 正常文库 `Nonimmune` 细胞（绝不触碰测试折患者），筛选 `|mean_NL - mean_NC| >= 0.5` 或 `max(mean_NL, mean_NC) >= 1.5` 的实质高表达/器官特异基因并与预设肝/结肠/分泌型 Ig 标志物求并集，在 top-512 HVG 筛选前将其方差置为 `-inf`，逐折记录被置换基因清单（`PCA/organ_marker_mask_audit.csv`）。
- **冻结 whole-human scGPT512 编码**：
  - 经 `canonical_to_target_vocab` 建立 `18,815` 核基因到 `scGPT vocab.json` 的严格单射映射（命中 `18,717` 个无碰撞词表 token），在共同核基因轴 `log1p(CP10k)` 上执行原版 51-bin 分位数分箱、`<cls>` 拼接、`max_length=1200` 确定性逐细胞随机截断与 L2 归一化。
  - 在 Intel Arc A750 XPU 上完成 `>= 64`（实测 `72`）个真实细胞的 CPU vs XPU 与跨批大小精度核验（同时记录 `min_cosine` 与 `max_abs_diff`）后编码全量 `70,628` 个肿瘤免疫细胞。

## 3. 池化对比、细类去混杂与 200 次患者置换零分布

- **池化与损失**：PCA32 与 scGPT512 分别比较均值池化（`mean`）与 `starter/linear_attention_mil.py` 透明线性 attention（`lambda_pair=0.1, lambda_l2=0.01`；附 `no_pair` 消融及 `SEED, SEED+1, SEED+2` 三种子归因稳定性）。
- **细类比例 vs 细胞内状态解耦**：
  - `PCA_subtype_balanced`：利用 M2 零 `unresolved` 分层细类，在每谱系两大主导细类（`T/NK`: `CD8 T cell` + `CD4 T cell`；`Mono/Mac/DC`: `Macrophage` + `Monocyte`；`B/Plasma`: `plasma cell` + `B cell`，均满足每类 `>= 5` 细胞、覆盖 100% 合格袋）内保留全部细胞求质心后严格 `50/50` 等权平均，彻底消除细胞饥饿与 PT/LM 亚型比例偏移。
  - `PCA_pure_primary` / `PCA_pure_secondary`：仅保留单一纯细类（`>= 20` 细胞）细胞建袋，直接测量纯单亚型转录状态的可分性。
- **200 次患者级符号翻转零分布**：在每折 `train_donors` 内执行 `200` 次独立患者级 `PT <-> LM` Rademacher 随机置换（`PCA/permutation_null_distribution.csv` & `permutation_null_summary.csv`），给出经验零分布均值、标准差、95% 分位区间与精确置换 $p$ 值。
