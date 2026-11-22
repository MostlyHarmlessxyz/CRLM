# CRLM 患者配对免疫程序与器官背景研究

围绕结直肠癌原发灶（PT）与肝转移灶（LM），本项目在严格患者级独立单位（`donor_uid`）与密码学预注册冻结下，解析跨患者、跨队列可复现的免疫细胞状态与分子程序，分离亚群组成变化与正常结肠–肝脏驻留背景，并基准评估线性表征与零样本单细胞基础模型（scGPT）及多实例学习（MIL）的外部泛化能力。原始研究需求见 [说明.md](docs/source_brief/说明.md) 与 [CRLM实验方案.md](docs/source_brief/CRLM实验方案.md)。

## 队列概览（质控后共 483,299 个体内单细胞 + 4 张空间切片）

| 数据集 | 研究角色 | 有效样本规模 | 质控保留细胞数 |
|---|---|---|---|
| `GSE245552` | 唯一主发现与模型开发集 | 16 位患者（13 对 PT–LM，39 个文库） | 173,450 |
| `GSE164522` | 外部复现与同患者四组织背景参照（PT/LM/NC/NL） | 10 位患者（40 份四组织样本） | 156,511 |
| `GSE178318` | 冻结后解锁的零表达暴露外部验证集（物理隔离 PBMC） | 6 位同步切除 IV 期 MSS 患者（12 个实体瘤文库） | 105,813 |
| `GSE203608` | 次要体内实体组织方向支持（物理隔离类器官） | 15 位体内患者（4 对 PT–LM） | 47,525 |
| `GSE335552` | 亚细胞原位空间微生态位验证（10x Xenium） | 3 位患者（3 张肝转移切片 + 1 张邻近肝切片） | 256,876 |

> 注：`GSE225857` 因公开元数据缺失可核实的患者级 PT–LM 配对映射，按预审计规则排除在主配对推断之外。

## 核心结论

- **方法学去偏与假阳性证伪**：构建官方 HGNC 单射基因桥（统一 18,815 / 18,733 共同核基因空间）与 20 桶检出率匹配去偏 UCell，消除肿瘤矩阵缺失 353 个线粒体/核糖体基因及秩截断测序深度偏倚；三队列一致证明初筛看似下调的 B/浆细胞未折叠蛋白反应（UPR）主要由浆细胞比例收缩（65.6%–77.3%）驱动，纯浆细胞配对伪批量无显著差异基因。
- **黏膜 Th17 程序逆转正常肝背景显著下调**：在同患者四组织队列 `GSE164522` 中，正常肝基线全员高于正常结肠（10/10 例），而肝转移灶逆转该梯度下调，通过器官背景家族校正（四联体观测交互 $I = -0.049$，Holm $P = 0.043$）；纯 CD4+ T 细胞在 `GSE164522`（9/10 例）与冻结后盲测集 `GSE178318`（5/6 例）一致同向下调。
- **脂质/瘢痕相关巨噬细胞重塑与跨队列模型基准**：核心脂质/瘢痕相关巨噬细胞模块（*SPP1*、*APOC1*、*APOE*、*GPNMB*、*TREM2*）在开发集（11/11 例）与盲测集 `GSE178318`（5/6 例）升高。15 个冻结分类头在 `GSE164522` 与 `GSE178318` 实现 100% PCA 特征覆盖与零拒判迁移：折内 32 维线性 PCA 在外部泛化与跨种子归因稳定性上全面匹配或优于 512 维零样本 scGPT，开发集上与 MixMIL 并列最优，HA/CTA 无额外增益。
- **Xenium 多尺度空间结构与功能验证方案**：在 `GSE335552` 切片中通过转录本优势比救回 3,584 个边界接触态 CD8+ T 细胞，`400 → 50 µm` 多尺度网格条件化表明粗尺度下 *SPP1*+ 髓系细胞与 T 细胞的距离差主要反映中尺度肿瘤巢分隔而非均匀近程排斥；据此制定分阶段人体原位杂交与自体髓系趋化验证方案（本研究已完成实际湿实验数为 0）。

## 快速入口

- **结果与复算**：[终稿总报告](artifacts/m8/final_report.md) · [复算指南](artifacts/m8/reproduce.md) · [局限与偏差说明](artifacts/m8/limitations_and_deviations.md) · [三版本对照表](artifacts/m6/three_version_evidence_comparison.md)
- **治理与状态**：[项目状态台账](STATUS.md) · [M5 密码学冻结记录](artifacts/m5/freeze_record.md) · [协作与防泄漏规范](AGENTS.md)
- **设计与准备**：[完整研究设计](docs/CRLM_Research_Plan.md) · [阶段实施计划](docs/IMPLEMENTATION_PLAN.md) · [算力指南](docs/COMPUTE_GUIDE.md) · [环境锁说明](environments/README.md)

## 仓库结构

| 目录 / 文件 | 说明 |
|---|---|
| `configs/project.json` | 预注册冻结锁定的队列分工、路径与超参数配置（`CRLM_V1.1_245_F1`） |
| `src/crlm/` | 核心 Python 库：HGNC 单射对齐、层级注释、去偏 UCell、多实例袋构建与先验 MIL |
| `scripts/` | `M0–M8` 全流程分析脚本、技术勘误脚本与自动化回归测试套件 |
| `starter/` | 经 106 项合成检查验证的配对效应、对称组成分解、Holm 校正与线性注意力原型 |
| `environments/` | 四个仓内 Python 3.12 虚拟环境（`.venv-m1/m2/m4/scgpt`）的精确 `uv` 锁文件 |
| `artifacts/m0/` … `m8/` | 各阶段正式统计表、阶段报告、出版级图表与 `output_sha256.txt` 清单 |
| `artifacts/archive/` · `techfix_v1_f1/` | 初版不可变快照（`CRLM_V1_245_F1`）与纯技术勘误分步消融对照轨 |
| `references/` | 85 条参考文献清单与锁定的上游对照源码（`MixMIL`、`hier-mil`） |
| `data/` · `embeddings/` · `models/` | 本地公开数据矩阵、A750 XPU `scGPT` 表示缓存与模型权重（大体量文件按 `.gitignore` 忽略） |
| `SHA256SUMS.txt` | 全仓受跟踪文件的 SHA-256 完整性校验清单 |

## 快速校验与回归测试

环境重建步骤见 [environments/README.md](environments/README.md)。在仓库根目录运行以下只读检查：

```bash
# 1. 结构完整性与冻结输入防篡改校验
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/check_project.py
sha256sum -c --quiet artifacts/m5/frozen_inputs_sha256.txt
sha256sum -c --quiet SHA256SUMS.txt

# 2. XPU 后端探针与全套方法学回归测试（含 starter 106 项合成自检）
PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python scripts/preflight.py --torch-probe
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/test_audit_remediations.py
```

> **运行硬件基线**：Intel Core i5-13400（64 GB RAM）+ Intel Arc A750 XPU（`torch 2.15.0+xpu`，相对 CPU 最大绝对误差 `2.09 × 10⁻⁷`）。公开数据、第三方代码与文献遵循各自原始许可条款。
