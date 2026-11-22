# CRLM 免疫程序项目

**研究目标：** 寻找结直肠癌原发灶（PT）与肝转移灶（LM）之间跨患者、跨独立来源可复现的免疫细胞状态与分子程序，分离细胞组成、亚型内转录状态与正常结肠–肝脏器官背景。项目原始需求以 [说明.md](docs/source_brief/说明.md) 为最高依据。

**[主交付报告与四幅主图（`CRLM_V1.1_AUDIT_REMEDIATED`）](artifacts/m8/final_report.md)**。本项目采用透明三版本并列架构：主线 **`CRLM_V1.1_AUDIT_REMEDIATED`**（冻结标识 `CRLM_V1.1_245_F1`）完整修复了单细胞排序深度偏倚、HGNC 基因别名脱节与粗谱系混杂；原版基线 **`CRLM_V1_245_F1`**（`M3–M8` 阶段产物与模型/嵌入快照）封存于 [archive/](artifacts/archive/CRLM_V1_245_F1/)；纯客观技术勘误 **`CRLM_V1_245_F1_techfix`** 对照存于 [techfix_v1_f1/](artifacts/techfix_v1_f1/)。核心实证结论：
- **`T/NK::B07_Th17_TWEAK` 显著负向四组织观测交互**：开发集 `GSE245552` `10/12` 配对患者下调（状态占比 `95.2%`）；主外部 `GSE164522` 同患者四组织观测交互 **`I = −0.04918`（`1+/9−`, `N=10`, `raw p = 0.0215`, 家族校正 `Holm p = 0.0430 <= 0.05` 通过）**，焦点纯 `CD4 T cell` 肿瘤差 `−0.02938`（`1+/9−`），显示结直肠原发肿瘤相对正常结肠正向偏移（`+0.03404`）而肝转移相对正常肝反向偏移（`−0.01514`），显著超出正常结肠–肝脏本底差；
- **`Mono/Mac/DC::Core_Lipid_SAM_R` 巨噬脂质/重塑核心程序**：开发集 `11/11` 患者升高（`+0.18092, p = 0.00098`），外部 `GSE164522`（`5+/3−`，四组织 `I` `4+/1−`；共线去污染 `Mac_SPP1_clean` `7+/1−` `D_T`、`5+/0−` `I`）与 `GSE203608` 纯巨噬（`3+/1−`，宽髓系 `2+/2−`）未达主族显著；初版 `B/Plasma::UPR` 经深度去偏、伪批量（`0` DEG）与对称分解（`72.6%–77.3%` 比例驱动）证伪淘汰；
- **`15/15` 零特征拒判全袋迁移与多尺度原位空间审计**：共享核基因空间下 `GSE164522` 全部 `15/15` 迁移分类头激活（6 个 PCA 头 `512/512 = 100%` HVG 覆盖），透明 `PCA` 在三谱系点估计上一致高于对应 `scGPT`；`GSE335552` 原位切片救回 `3,584` 个接触态 `CD8 T`（其中 `M397T` `3,236` 个）后显示 `200–400 µm` 距离差随尺度细化大幅衰减，且在 `M397T` 主要集中于距 `EPCAM+` 上皮 `<15 µm` 肿瘤巢带。[三版本对照总表](artifacts/m6/three_version_evidence_comparison.md) · [复算指南](artifacts/m8/reproduce.md) · [项目状态](STATUS.md)

## 从这里开始

1. 研究问题、统计框架与阶段门禁见 [完整研究设计](docs/CRLM_Research_Plan.md) 与 [实施计划](docs/IMPLEMENTATION_PLAN.md)。
2. 患者身份核实、本机 106 项基线自检与 `scGPT` 预训练重叠审计见 [M0 报告](artifacts/m0/stage_report.md)；全开发集 39 库（`173,450` 保留细胞，零 `unresolved` 层级注释）见 [M2 报告](artifacts/m2/stage_report.md)。
3. 去偏 `UCell` 候选筛选与配对伪批量见 [M3 报告](artifacts/m3/stage_report.md)；全袋分类消融与 A750 XPU 精度验证见 [M4 报告](artifacts/m4/stage_report.md)；冻结契约与历史暴露披露见 [M5 报告](artifacts/m5/stage_report.md)；外部队列与空间/实验设计见 [M6 报告](artifacts/m6/stage_report.md) 与 [M7 报告](artifacts/m7/stage_report.md)。

## 目录与用途

| 路径 | 用途及当前状态 |
|---|---|
| `AGENTS.md` | 持续生效的科学、统计与防泄漏工程约束 |
| `docs/` | 原始需求（`source_brief/说明.md`）、68 页研究计划、M0–M8 实施计划、文献准备与算力指南 |
| `configs/project.json` | 与 M5 冻结快照保持同哈希锁定的项目配置（`V1.1_AUDIT_REMEDIATED` / `CRLM_V1.1_245_F1`） |
| `src/crlm/` | 核心可复用库：官方 HGNC 单射对齐、两阶段层级注释、20 桶检测率匹配去偏 `UCell`、全袋构建与先验 MIL |
| `starter/` | 经 106 项合成检查验证的配对效应、对称组成分解、Holm 注册校正与透明线性 attention 原型 |
| `scripts/` | `M0–M8` 全链路实证分析、Track 1 纯技术勘误与自动化回归测试脚本（`test_audit_remediations.py`） |
| `references/` | 85 条主参考文献元数据、本地论文索引及仓内锁定的 `MixMIL` / `hier-mil` 第三方源码（`references/code/`） |
| `data/`、`embeddings/`、`models/` | 本地真实单细胞对象、A750 XPU `scGPT` 表示缓存、全集迁移模型及 `CRLM_V1_245_F1` 历史归档 |
| `artifacts/m0/` … `artifacts/m8/` | `CRLM_V1.1_AUDIT_REMEDIATED` 各阶段正式报告、结果表、主图与 `output_sha256.txt` 清单 |
| `artifacts/archive/CRLM_V1_245_F1/` | 初版 `M3–M8` 阶段产物与根级快照及 `archive_sha256.txt` |
| `artifacts/techfix_v1_f1/` | Track 1 冻结程序与冻结分类器纯技术勘误分步消融结果 |

## 快速校验与全套回归测试

```bash
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/check_project.py
PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python scripts/preflight.py --torch-probe
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/test_audit_remediations.py
```

## 当前硬件实测结论

本机 Intel Core i5-13400（64 GB RAM）与 Intel Arc A750 XPU（隔离 `.venv-scgpt`，`torch 2.15.0+xpu`，72 细胞实测相对 CPU 最大绝对误差 `2.09 × 10⁻⁷`）已完成开发集 `39` 个文库（`173,450` 保留细胞，其中肿瘤全袋 `70,628` 个免疫细胞）与外部队列 `GSE164522`（`156,511` 洁净细胞，其中肿瘤全袋 `68,050` 个免疫细胞）、`GSE203608`（`47,525` 洁净体内细胞）及 `GSE335552`（`4` 张原位切片）的全量计算与哈希闭环验证。`GSE178318` 保持零表达矩阵暴露（仅 M0 核查公开 GEO 临床/配对元数据），待获取作者书面使用条款；论文、代码与公开数据的原有许可条款持续适用。
