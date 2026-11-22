# 方法学局限、审计整改记录与方案偏差说明（`CRLM_V1.1_AUDIT_REMEDIATED`）

## 1. 相对研究计划（`docs/CRLM_Research_Plan.md`）的受控工程实现调整

1. **配对伪批量负二项检验实现（Python 替代 R `edgeR`）**：
   - 原计划第 7.2 节建议优先采用 `edgeR::glmQLFTest` 准似然 $F$ 检验；鉴于目标主机本地工作流统一运行于隔离 Python 环境，本研究在 `scripts/analyze_paired_counts.py` 中采用 `PyDESeq2` 拟合 `~ donor_uid + tissue` 配对负二项 GLM（显式关闭 `N=2/donor` 配对设计下失真的 Cook's 距离过滤，对 Wald/$t$ 结果严格屏蔽 `_genewise_converged / _MAP_converged / _LFC_converged` 任一未收敛基因），并将默认渐近正态 Wald 检验升级为**有限样本残差自由度 `df = N_pairs - 1` 双侧 Student-$t$ 检验**；同时并列在供者内可交换性假设下对 `log₂(CP10k + 0.5)` 配对差精确遍历半空间 `32–2,048` 个符号向量（对应对称全空间 $2^{N_{\text{pairs}}} = 64\text{–}4,096$ 种翻转配置）计算双侧精确符号翻转检验（`pvalue_signflip`）。在无信号谱系 `B/Plasma`（`df=7`）与纯 `plasma cell`（`df=5`）上，Student-$t$ 与符号翻转均给出 `0` 个 `BH < 0.05` 基因（对比同一 `V1.1` 输入下渐近正态 Wald 的 `42` 与 `39` 个显著呼叫及原版归档 `B/Plasma` 的 `47` 个呼叫），在髓系（`661` vs `611`）与纯巨噬（`400` vs `399`）上高度一致。[配对计数汇总](../m3/paired_counts/fit_summary.csv)
2. **Intel Arc A750 XPU 硬件后端与确定性边界**：
   - 原计划针对无 CUDA 环境以 CPU 为基线、XPU 待实测验证。本机采用隔离 `.venv-scgpt`（**`torch 2.15.0+xpu`**）在 `72` 个真实三谱系肿瘤细胞上完成与 CPU `float32` 的逐细胞对位核验：最小余弦相似度 **`0.99999982`**、最大逐维绝对误差 **`2.09 × 10⁻⁷`**，并在 `16` 个细胞上完成逆序跨批大小重排核验（最大绝对误差 **`8.94 × 10⁻⁸`**，见 [scgpt_profile.json](../m4/resource_profiles/scgpt_profile.json)），在此实证基础上完成开发集 `70,628` 个与外部 `68,050` 个细胞的快速编码。超过 `1,199` 表达基因的细胞通过基于 `cell_uid` SHA256 派生的独立随机种子抽样，完全不依赖批次或读入顺序。

## 2. 三版本演进中已查明并闭环修复的核心缺陷记录

1. **单细胞程序评分深度偏倚修复**：初版自定义 `positive_midrank_1500` 仅对正表达基因取秩并令零表达为 `0`，在稀疏 scRNA-seq（中位检出基因 `< 1,500`）下退化为检测深度的递增函数（合成零效应 Poisson 检验 `max |ρ| = 0.892`）。`V1.1` 升级为全测量轴 Mann–Whitney `UCell` 减去开发集锁定的 **20 桶 × 10 重复检测率匹配背景对照**（`debiased_ucell_scores`，合成检验残余相关 `max |ρ| = 0.0755 < 0.08`），并在外部队列严格复用冻结 `matched_controls.csv`，杜绝外部自适应重拟合。
2. **跨队列 HGNC 基因轴与共享核空间对齐**：通过官方 HGNC 快照（含基于 `ENSG` 的重名符号 `SEPT7/SEPTIN7` 单射消歧）构建 `245_164_shared_nuclear`（`18,815` 核基因）与 `245_164_203608_shared_nuclear`（`18,733` 核基因），显式剥离 `GSE245552` 仅在邻近组织出现的 `353` 个特征（含全部 13 个 `MT-` 基因及若干核糖体蛋白、激酶与非编码转录本），在 Track 1 原版冻结程序上恢复了 `1–8` 个因别名丢失的成员（`UPR` `+1`、`Mac_SPP1` `+6~+8`）并使 `T/NK`、`Mono/Mac/DC`、`B/Plasma` 全部 6 个外部 PCA 模型达到 **`512/512 = 100.0%`** 特征覆盖。
3. **两阶段无冲突层级注释与跨队列细胞清洗**：解耦 `CSF3R + S100A8/A9`（`FCGR3B=0`）单核/巨噬误判、补齐 `CD4/CD8` 层级拆分使 `M2` 已解析谱系 `unresolved` 率降至 `0.00%`；同步剔除 `GSE164522` 全部 `6,836` 个跨谱系 `CD3E+` 双细胞群（含初版漏检的 `hC75_B-CD3E` `3,174` 个细胞）与 `GSE203608` 的 `12,136` 个 `Cell_type=='NA'` / 高线粒体死细胞，并在两外部队列提供 `annotation_m1.classify` 同规则重注释对照分支。
4. **先验 MIL 实现与评估去偏**：本地 vendor 锁定 `MixMIL` 与 `hier-mil` 源码，修复 `MixedMIL` 初始化投影缩放奇偶切片 `[::2]` 丢弃半数投影值与漏除 `1/sqrt(Q)` 导致的先验方差膨胀，将 `HierarchicalMIL` 升级为严格折内 3 折供者分层内层 CV 调参，并以全袋评估 + 200 次患者置换零分布替代初版单次 256 细胞截断与仅翻转 3/13 对患者的偏倚哈希压力测试。

## 3. 固有科学边界与已量化的技术边界

- **既往外部暴露与确认属性边界**：因初版 `CRLM_V1_245_F1` 已访问 `GSE164522` 与 `GSE203608` 表达矩阵，且共享核基因轴 `245_164_shared_nuclear` 使用了外部公开基因符号列表确定可迁移特征范围，本研究严格遵守 `AGENTS.md` 第 8 节，将 `V1.1` 在 `164/203608` 上的检验界定为**已披露历史暴露的外部复现与稳健性重评估**；真正的首次未暴露表达矩阵外部检验留待获取 `GSE178318`（`N=6` 对，M0 仅核查公开元数据）PI 书面使用条件或新队列后执行。
- **样本量与离散检验功效边界**：独立推断单位严格为真实患者 `donor_uid`（`GSE164522` 配对 `N = 8–10`、髓系四组织 `Q` 子集 `N = 5`；`GSE203608` 体内配对 `N = 2–4`，且 `Patient 49` 临床分期注释存在文献内不一致，已单列剔除敏感性；`GSE225857` 因缺供者映射继续排除）。在 `N = 5` 时双侧精确符号检验即使 `5/0` 全同向理论最小 `raw p = 0.0625 > 0.05`，且 `GSE164522` 四组织中 `Core_Lipid_SAM_R`（`I = +0.07425, 4/1`）与 `Resident_Kupffer`（`I = +0.06798, 4/1`）均呈小样本正向，驻留巨噬背景尚不能单凭现有 `N=5` 完全排除。
- **`scGPT` 预训练与表示解释边界**：语料重叠审计（[scgpt_pretraining_overlap_audit.json](../m0/scgpt_pretraining_overlap_audit.json)）将 `GSE245552` 评定为 `DEFINITE_ZERO_OVERLAP`、`GSE203608` 评定为 `NEGLIGIBLE_POST_CUTOFF_PUBLIC_RELEASE`、`GSE164522` 评定为 `POTENTIAL_UNSUPERVISED_EXPOSURE_UNVERIFIED`（公开权重不含逐供者哈希清单）；高组织分类 AUROC 反映位点微环境转录可分性而非转移因果驱动。
- **已独立量化的两项细微边界（红队复核存档）**：
  1. `program_scores.py` 对照抽样在极小候选池桶内的同重复回退：经全表扫描 `artifacts/m3/program_definitions_draft/matched_controls.csv`，**全部 19 个 M3 筛选集（18 个候选程序 + 1 个器官驻留对照，含两个冻结主候选 `Core_Lipid_SAM_R`、`B07_Th17_TWEAK` 及 `B19_S4_Mac_SPP1_clean`）同重复内重复对照数均为 `0`**，仅非肿瘤背景参照大集 `B06_S6_SAM`（`651` 个被测基因 / `6,510` 次抽样）含 `150` 次同重复回退；
  2. `annotation_m1.py` 细类与面板敏感性：阶段二 `APOE` 仅影响开发集 `1.4%–1.7%` 巨噬细类划分（不依赖 `APOE` 的粗类 `Mono/Mac/DC` 与 `monocyte/macrophage` 在开发集均为 `11/11 = 100%`）；在 `GSE164522` 处理后连续表达矩阵上，纯 `Macrophage` 的 `Core_Lipid_SAM_R` 在原作者聚类标签下为 `+0.02847 (5+/3−)`，而在简化阈值规则 `rule_annotation_m1` 下为 `−0.00589 (3+/5−)`，完整并列披露于 [focal_subtype_effects_all_cohorts.csv](../m6/focal_subtype_effects_all_cohorts.csv)；`ASGR1+ASGR2` 双阳在 M2 共见于 `544 / 173,450 (0.31%)` 个保留细胞（若关闭双阳肝面板贡献仅改变 `423 / 173,450 = 0.24%` 个粗谱系标签）。
- **空间面板覆盖与观测推断边界**：`I = D_T_Q − D_N_Q` 是四组织观测交互项而非因果转移效应；`GSE335552` 仅含肝转移与邻近肝切片（无配对原发结直肠切片及配对 H&E 缺氧/坏死边界），空间检验为切片内局部方格置换，机制判别须依 [M7 组织原位与同供者趋化验证预案](../m7/candidate_validation_design.md) 执行。
