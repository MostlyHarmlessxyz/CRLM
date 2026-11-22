# `scGPT` 预训练语料与 CRLM 项目队列重叠度专项核查报告

> **合规依据**：`AGENTS.md` 第 7 节（“scGPT 缓存记录权重/词表/代码/基因背景/输入变换/逐细胞随机种子；冻结权重不保证无预训练暴露”）及 `docs/IMPLEMENTATION_PLAN.md` M4 验收规范。

## 1. 本地模型权重、词表与源码溯源

本项目在 M4（开发集细胞表示比较）与 M6（`GSE164522` 冻结肿瘤分类器外部评估）中仅以**冻结推理模式（零梯度更新、零微调）**调用 `scGPT-whole-human` 检查点：

| 资产文件 | 路径 | 字节数 | SHA256 | 关键元数据 |
|---|---|---:|---|---|
| 模型权重 | `models/scGPT-whole-human/best_model.pt` | `205,385,258` | `6cb5d451ab5c4b33eb673adbe4fddc61d2389df1b89b7651a9fe2e557572b922` | HuggingFace `wanglab/scGPT-human` (`rev a24c2377...`)，12 层 Transformer，`embsize=512`，`n_bins=51` |
| 词表文件 | `models/scGPT-whole-human/vocab.json` | `1,317,639` | `acca93d114ca62c3f0f50debbd23e8c87f0714f4737764454f6b2b13f2e8580f` | 共 `60,697` 个 token（含 `<pad>`、`<cls>`、`<eoc>`） |
| 训练参数 | `models/scGPT-whole-human/args.json` | `1,300` | `c18e075e018140cb8b2d9029387b9de26607a5ce6a8ccabd6ead70cd76b95d60` | `save_dir`: `.../cellxgene_census_human-May23-08-36-2023` |
| 上游源码 | `models/scGPT-whole-human/scGPT_source.tar.gz` | `12,072,234` | `46be386b0fae377844f908a007cdda2a7c29c75667e8460b88879b60a677bb81` | 上游提交 `cebd6fae655b9c585a4807daa3ac31bb764f06b4` |

## 2. 预训练语料时间截点、器官/疾病筛选与监督信号核查

核查 `models/scGPT-whole-human/args.json` 与上游源码 `models/scGPT-whole-human/scGPT-cebd6fae655b9c585a4807daa3ac31bb764f06b4/data/cellxgene/{data_config.py,cancer_list.txt,query_list.txt,metainfo.json,README.md}` 可确证以下四项事实：

1. **CELLxGENE Census 快照硬时间截点（`2023-05-08`）**：
   - `data_config.py` 第 5 行明确定义 `VERSION = "2023-05-08"`；
   - `args.json` 第 2–3 行记录 `"data_source": "/scratch/ssd004/datasets/cellxgene/scb_strict/human"` 与 `"save_dir": "/scratch/ssd004/datasets/cellxgene/save/cellxgene_census_human-May23-08-36-2023"`。
   - 因此，任何在 **2023 年 5 月 8 日之后**才公开释放的单细胞数据集，物理上不可能进入该检查点的预训练语料。
2. **`whole-human` 与 `pan-cancer` 筛选查询范围**：
   - 本地检查点对应 `scb_strict/human`（即 `whole-human` 正常细胞主语料）：`data_config.py` 第 4、12–19 行对七大器官（`MAJOR_TISSUE_LIST = ["heart", "blood", "brain", "lung", "kidney", "intestine", "pancreas"]`）及 `"others"` 分区均强制施加 **`suspension_type != 'na' and disease == 'normal'`**。其中结肠/肠道正常细胞命中 `tissue_general == 'intestine'`，外周血命中 `'blood'`，正常肝组织（`'liver'`）命中 `'others'` 分区。
   - 同目录 `cancer_list.txt` 第 17 行显式包含 **`colorectal cancer`**（第 15 行为 `adenocarcinoma`），用于 `query_list.txt` 第 9 行的 `pan-cancer` 查询分支（`disease == 'colorectal cancer'` 等 21 类肿瘤）。
   - 同目录 `metainfo.json` 列出的 60 个早期 CELLxGENE Portal 集合中，经全文检索与 `GSE245552`、`GSE164522`、`GSE178318`、`GSE203608`、`colorectal`、`liver metastasis` 均**零直接命中**；但由于 CELLxGENE Census `2023-05-08` TileDB-SOMA 全库未随模型权重附带完整 `dataset_id` 清单，早于 `2023-05-08` 公开且被社区收录进 CELLxGENE 的研究仍可能通过邻近正常组织细胞（标注 `disease == 'normal'`）或跨检查点预训练管线产生无监督表示接触。
3. **预训练任务无监督属性（零 `PT` vs `LM` 标签泄漏）**：
   - `args.json` 第 31–32、45–48 行明确记录 `"no_cls": true`、`"USE_CLS": false`、`"no_cce": true`、`"USE_CCE": false`、`"MVC": true`、`"USE_GENERATIVE_TRAINING": true`。
   - 预训练目标纯属单细胞内部分位数分箱表达值的掩码重建（mask ratio `0.25 / 0.50 / 0.75`），**完全不含患者 ID、肿瘤原发灶 vs 肝转移灶（`PT` vs `LM`）分类标签或配对对比信号**。

## 3. 项目各队列重叠度逐项判定

| 数据集 | 项目角色 | GEO 提交日期 | GEO 公开日期 / 论文发表 | 相对 `2023-05-08` 截点 | 筛选规则潜在命中分析 | 重叠风险评级 | 对本研究推断的约束与解释边界 |
|---|---|---|---|---|---|---|---|
| **`GSE245552`** | V1/V1.1 唯一开发集（M1–M5） | `2023-10-17` | `2024-01-05` / *Cancer Res* 2024 `[B07]` | **严格晚于截点**（连首次提交也晚 5 个月以上） | 物理上不可能存在于 `2023-05-08` Census | **`DEFINITE_ZERO_OVERLAP`** | M4 开发集折外（OOF）`scGPT` vs `PCA` 比较完全不存在预训练语料重叠 |
| **`GSE164522`** | 主外部复现与四组织背景（M6） | `2021-01-10` | `2022-02-12` / *Cancer Cell* 2022 `[B02]` | **早于截点**（公开早约 15 个月） | 若在 `2023-05-08` 前被 CELLxGENE Census 收录：配对正常结肠 `PN`（`intestine`）、正常肝 `MN`（`others: liver`）与 `PBMC`（`blood`）命中 `disease == 'normal'`；肿瘤细胞在 `pan-cancer` 配置中命中 `colorectal cancer` | **`POTENTIAL_UNSUPERVISED_EXPOSURE_UNVERIFIED`** | ① 主生物学复现检验（`Mac_SPP1`、`UPR` 等患者级配对效应）完全基于非模型基因集评分，**不依赖 `scGPT`**；② M6 模型迁移必须与零预训练的同折同细胞 **`PCA` 基线**并列报告，不得将 `scGPT` AUROC 单独夸大为“零先验暴露的纯盲外部泛化” |
| **`GSE178318`** | 次要外部队列（当前 `terms_pending` 未下载表达） | `2021-06-16` | `2021-09-10` / *Cell Discov* 2021 `[B03]` | **早于截点**（公开早约 20 个月） | 含原发 CRC、肝转移与 PBMC；若被 Census 收录，正常标注细胞或肿瘤细胞存在与 `GSE164522` 同类的潜在无监督收录可能 | **`POTENTIAL_UNSUPERVISED_EXPOSURE_UNVERIFIED`（当前未启用）** | 当前維持 `terms_pending` 未运行任何编码；未来若获作者确认接入，主验证同样仅用非模型基因程序统计，若跑 `scGPT` 须附同样声明与 `PCA` 对照 |
| **`GSE203608`** | 可选次要外部队列（M6） | `2022-05-23` | `2023-10-20` / *Nat Commun* 2023-11 `[B46]` | 提交早于截点，但 **GEO 公开晚于截点 5.4 个月** | `2023-05-08` 时仍处未公开保密状态，极不可能进入公开 CELLxGENE Census；且 `metainfo.json` 零命中 | **`NEGLIGIBLE_POST_CUTOFF_PUBLIC_RELEASE`** | 本研究对 `GSE203608` 仅计算基因程序配对方向与幅度，不运行 `scGPT` 模型推理，无表示泄漏影响 |
| **`GSE225857`** | V1/V1.1 排除（患者身份映射待核） | `2023-02-22` | `2023-02-22` / 2023 | 早于截点 2.5 个月 | 潜在可被 Census 收录，但本阶段不纳入配对推断或模型评估 | **`NOT_USED_EXCLUDED`** | 维持排除状态 |
| **`GSE298084`** | 癌旁肝背景参考 | `2025-05-27` | `2025-08-11` / 2025 | **严格晚于截点**（晚 >2 年） | 物理上不可能重叠 | **`DEFINITE_ZERO_OVERLAP`（未用于 `scGPT`）** | 无重叠风险 |
| **`GSE335552`** | M7 10x Xenium 空间定位支持 | `2026-06-16` | `2026-08-03` / 2026 | **严格晚于截点**（晚 >3 年，且为 377 基因靶向空间原位成像） | 物理上不可能重叠 | **`DEFINITE_ZERO_OVERLAP`（未用于 `scGPT`）** | 仅用于 M7 空间邻近几何检验，不涉及 `scGPT` |

## 4. 工程与方法学落地要求

1. **公平基线绑定**：在 M4（`GSE245552`）和 M6（`GSE164522`）中，所有 `scGPT` 表示模型（`scGPT_mean`、`scGPT_attention`、`scGPT32`/`scGPT512` `HA/CTA/MixMIL`）必须与仅在 `GSE245552` 训练折拟合的 **`PCA` 表示模型**使用完全相同的患者划分、细胞集合与共同核基因背景（`SHARED_NUCLEAR_UNIVERSE`）。
2. **生物结论脱钩**：跨队列生物学复现结论完全由冻结候选程序的患者级配对检验（双侧精确符号检验 + Holm 校正、四组织交互项 `I = D_T - D_N`）决定；`scGPT` 分类器仅作为补充性的组织来源可分性与表示比较证据。
3. **透明披露**：在 M4、M6 与 M8 总报告中引用本核查文件（`artifacts/m0/scgpt_pretraining_overlap_audit.md` 及同名 `.json`），明确区分 `GSE245552` 的确定性零重叠与 `GSE164522`（2022）无法排除 CELLxGENE Census `2023-05-08` 无监督正常细胞收录的客观边界。
