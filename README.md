# CRLM 免疫程序项目

**研究目标：** 寻找结直肠癌原发灶（PT）与肝转移灶（LM）之间跨患者、跨独立来源可复现的免疫细胞状态与分子程序，分离细胞组成、亚型内转录状态与正常结肠–肝脏器官背景。项目原始需求以 [说明.md](docs/source_brief/说明.md) 为最高依据。

**[主交付报告与四幅主图（`CRLM_V1.1_AUDIT_REMEDIATED` · 8/8 预登记检验全观测闭环）](artifacts/m8/final_report.md)**。本项目采用透明三版本并列架构：主线 **`CRLM_V1.1_AUDIT_REMEDIATED`**（冻结标识 `CRLM_V1.1_245_F1`）完整修复了单细胞排序深度偏倚、HGNC 基因别名脱节与粗谱系混杂；原版基线 **`CRLM_V1_245_F1`** 不可变封存于 [archive/](artifacts/archive/CRLM_V1_245_F1/)；纯客观技术勘误 **`CRLM_V1_245_F1_techfix`** 对照存于 [techfix_v1_f1/](artifacts/techfix_v1_f1/)。核心实证结论：
- **`T/NK::B07_Th17_TWEAK` 显著负向四组织观测交互与前瞻复现**：开发集 `GSE245552` `10/12` 配对患者下调（状态占比 `95.2%`）；主外部 `GSE164522` 同患者四组织观测交互 **`I = D_T_Q − D_N_Q = −0.04918`（`1+/9−`, `N=10`, `raw p = 0.0215`, 背景族 `Holm p = 0.0430 <= 0.05` 通过）**，焦点纯 `CD4 T cell` 在 `164`（`−0.02938, 1+/9−`）及 M5 冻结后首次解锁表达矩阵的**前瞻零表达暴露队列 `GSE178318`（`−0.02067` `[−0.04040, −0.00169]`, `1+/5−`；宽谱系 `−0.00986, 2+/4−`）**中均高比例同向下调；
- **`Mono/Mac/DC::Core_Lipid_SAM_R` 巨噬脂质/重塑核心程序**：开发集 `11/11` 患者升高（`+0.18092, p = 0.00098`），**在零表达暴露前瞻队列 `GSE178318` 达 `5+/1−` 升高（`+0.09936`，纯巨噬 `+0.10183` `[+0.02595, +0.19140]`, `5/1`）**，在 `GSE164522`（`5+/3−`，四组织 `I` `4+/1−`；共线 `Mac_SPP1_clean` `7+/1−` `D_T`、`5+/0−` `I`）与 `GSE203608` 纯巨噬（`3+/1−`）中全部同向；初版 `B/Plasma::UPR` 跨独立队列一致证伪淘汰；
- **双外部队列（`164` + 前瞻 `178`）各 `15/15` 零拒判迁移与多尺度空间微解剖审计**：两个外部队列均实现 `512/512 = 100.0%` PCA HVG 覆盖，零样本 `scGPT` 未相对简约线性 `PCA` 提供一致外部增量（`164` 三谱系 `PCA` 均不劣于或优于 `scGPT_mean`；前瞻 `178` `T/NK`: `PCA 0.917` vs `scGPT_mean 0.583`, `Δ = +0.333 [+0.083, +0.556]`，髓系 `PCA_attention = 0.944, 6/6` 对，而全队列 `10/12` 对 CI 跨 `0` 且 `B/Plasma` 注意力 `scGPT 0.667 > PCA 0.625`）；`GSE335552` 在 `M397T` 救回 `3,236` 个（四片合计 `3,584` 个）接触态 `CD8 T` 后结合多尺度（`400->50 µm`）、半格平移与 `<15 µm EPCAM+` 分层量化了微解剖构型效应（湿实验执行数为 `0`）。[三版本对照总表](artifacts/m6/three_version_evidence_comparison.md) · [复算指南](artifacts/m8/reproduce.md) · [项目状态](STATUS.md)

## 从这里开始

1. 研究问题、统计框架与阶段门禁见 [完整研究设计](docs/CRLM_Research_Plan.md) 与 [实施计划](docs/IMPLEMENTATION_PLAN.md)。
2. 患者身份核实、本机 106 项基线自检与 `scGPT` 语料重叠审计见 [M0 报告](artifacts/m0/stage_report.md)；全开发集 39 库（`173,450` 保留细胞，已解析谱系零 `unresolved` 层级注释）见 [M2 报告](artifacts/m2/stage_report.md)。
3. 去偏 `UCell` 筛选与配对伪批量见 [M3 报告](artifacts/m3/stage_report.md)；全袋分类消融见 [M4 报告](artifacts/m4/stage_report.md)；预解锁冻结契约见 [M5 报告](artifacts/m5/stage_report.md)；**3 个外部队列共 8/8 项已测预登记检验（1 项四组织交互 `I` 通过 Holm 0.05，0/2 项主肿瘤 `D_T` 达显著；含 [GSE178318 前瞻零表达暴露报告](artifacts/m6/GSE178318/stage_report.md)）**与原位空间审计/验证方案见 [M6 报告](artifacts/m6/stage_report.md) 与 [M7 报告](artifacts/m7/stage_report.md)。

## 目录与用途

| 路径 | 用途及当前状态 |
|---|---|
| `AGENTS.md` | 持续生效的科学、统计与防泄漏工程约束 |
| `docs/` | 原始需求（`source_brief/说明.md`）、68 页研究计划、M0–M8 实施计划、文献准备与算力指南 |
| `configs/project.json` | 与 M5 预解锁冻结包同哈希锁定的项目配置（`V1.1_AUDIT_REMEDIATED` / `CRLM_V1.1_245_F1`） |
| `src/crlm/` | 核心可复用库（M5 锁定后 0 字节改动）：官方 HGNC 单射对齐、两阶段层级注释、去偏 `UCell`、全袋构建与先验 MIL |
| `starter/` | 经 106 项合成检查验证的配对效应、对称组成分解、Holm 注册校正与透明线性 attention 原型 |
| `scripts/` | `M0–M8` 全链路分析（含 `GSE178318` 前瞻验证）、Track 1 技术勘误与 7 组自动化回归套件（`test_audit_remediations.py`） |
| `references/` | 85 条主参考文献元数据、本地论文索引及仓内锁定的 `MixMIL` / `hier-mil` 第三方源码（`references/code/`） |
| `data/`、`embeddings/`、`models/` | 四个单细胞队列稀疏对象（3 份 `178` PBMC 物理隔离）、A750 XPU `scGPT` 表示缓存、全集迁移模型及 `V1` 归档 |
| `artifacts/m0/` … `artifacts/m8/` | `CRLM_V1.1_AUDIT_REMEDIATED` 各阶段正式报告、结果表、主图与 `output_sha256.txt` 清单 |
| `artifacts/archive/CRLM_V1_245_F1/` | 初版交付完整不可变快照与 `archive_sha256.txt` |
| `artifacts/techfix_v1_f1/` | Track 1 冻结程序（8/8 检验）与冻结分类器纯技术勘误分步消融结果 |

## 快速只读校验与全套回归测试

```bash
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/check_project.py
PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python scripts/preflight.py --torch-probe
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/test_audit_remediations.py
```

## 当前硬件实测结论

本机 Intel Core i5-13400（64 GB RAM）与 Intel Arc A750 XPU（隔离 `.venv-scgpt`，**`torch 2.15.0+xpu`**，72 细胞实测相对 CPU 最大绝对误差 `2.09 × 10⁻⁷`）已完成开发集 `GSE245552`（`173,450` 保留细胞，肿瘤全袋 `70,628` 细胞）、主外部 `GSE164522`（`156,511` 洁净细胞，肿瘤全袋 `68,050` 细胞）、前瞻外部 `GSE178318`（12 个实体瘤文库 `105,813` 保留细胞，编码肿瘤三谱系细胞 `59,338` 个、34 个 `>=50` 细胞合格袋入袋 `59,298` 个）、次要体内 `GSE203608`（`47,525` 洁净细胞）及原位空间 `GSE335552`（`4` 张 Xenium 切片）的全量计算与 SHA256 闭环验证。论文、代码与公开数据的原有许可条款持续适用。
