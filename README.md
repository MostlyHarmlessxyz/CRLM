# CRLM 免疫程序项目

**研究目标：** 寻找结直肠癌原发灶与肝转移灶之间跨患者、跨独立来源可复现的免疫细胞状态与分子程序，评估组成、状态及器官背景。项目意图以 [说明.md](docs/source_brief/说明.md) 为最高依据。

GSE245552 的全部 39 个公开样本已从原始交付计数导入稀疏 H5AD，覆盖 16 位患者。190,190 个输入细胞中有 190,031 个通过开发集 QC。患者与各谱系实际覆盖见 [M2 阶段报告](artifacts/m2/stage_report.md)，当前分析进度见 [STATUS.md](STATUS.md)。

## 从这里开始

1. 研究问题、数据角色和阶段验收见 [完整报告](docs/CRLM_Research_Plan.md) 与 [实施计划](docs/IMPLEMENTATION_PLAN.md)。
2. 患者身份和本机自检见 [M0 报告](artifacts/m0/stage_report.md)。重建 Pt18 对象、检查基因范围和初级谱系覆盖见 [M1 报告](artifacts/m1/stage_report.md)。
3. 患者程序分析使用 [M3 开发规则](artifacts/m3/program_definitions_draft/development_rules.md)，模型比较使用 [M4 患者分折方案](artifacts/m4/benchmark_design.md)。

## 目录与用途

| 路径 | 用途及当前状态 |
|---|---|
| `AGENTS.md` | 简明、持续生效的研究与工程约束 |
| `docs/source_brief/` | 用户原始材料，保留原文件名 |
| `docs/CRLM_Research_Plan.md`、`docs/report_assets/` | 对应原 68 页报告的正文、附录、参考文献和三张图 |
| `docs/IMPLEMENTATION_PLAN.md` | M0–M8 的输入、产物、验收及停止规则 |
| `docs/PAPERS_AND_DATA_PREPARATION.md` | 哪些论文先读、缺哪些 Methods/补充表、如何交给 Agent |
| `docs/COMPUTE_GUIDE.md` | i5-13400/64 GB/A750 本地路线与按需 GPU 判断 |
| `prompts/` | Plan 与 M0–M1 执行提示词，按使用阶段选择 |
| `configs/project.json` | **草案**配置；患者单位、数据角色、资源默认，不是冻结注册或已完成结果 |
| `starter/` | 原始、已验证的小型实现；直接看自己的 `README.txt`，不要当成完整单细胞流水线 |
| `scripts/preflight.py` | 无第三方依赖的只读环境/资源探测；可选保存本次结果 |
| `scripts/check_project.py` | 只读检查交接文件、相对路径、角色和参考表是否齐全 |
| `scripts/import_gse245552.py`、`src/crlm/annotation_m1.py` | GSE245552 的真实稀疏导入、重读核对与初级谱系标记 |
| `workflow/` | 阶段运行与输出约定；没有冒充可运行全流程的 DAG |
| `references/reading_manifest.csv` | 85 条文献记录的优先级、获取/阅读状态及用途 |
| `references/papers/` | 已归档的本地论文 PDF，见 [本地文件清单](references/local_files.csv) |
| `references/supplements/`、`references/code/` | 新增原始补充材料与第三方代码 provenance 的位置 |
| `data/raw/`、`data/processed/`、`data/manifests/` | 真实数据的位置；交付时没有新增全量矩阵 |
| `embeddings/`、`models/` | 冻结表示和模型版本的位置；交付时为空 |
| `artifacts/m0/` … `artifacts/m8/` | 各阶段的正式执行与验收记录 |
| `artifacts/handoff/` | **本次交接包核验**，不属于真实研究或用户主机验收 |

## 两条现在可运行的轻量命令

这些命令只读并向终端输出 JSON，不安装软件、不联网、不下载数据、不训练。需 Python 3.10+；原 `starter/` 的已验证环境为 Python 3.12.14。

```bash
python scripts/check_project.py
python scripts/preflight.py
```

脚本按自身位置寻找项目，换工作目录也能使用其完整路径调用。`preflight.py` 的运行位置就是被检测机器，**交付包内的沙箱检测不能代表你的本地主机**。

在执行阶段，若已安装 PyTorch 且希望显式探测其后端：

```bash
python scripts/preflight.py --torch-probe
```

本次主机检测保存在 [M0 环境记录](artifacts/m0/environment_report.md)。

CPU 基础环境、单细胞分析环境、scGPT 环境应隔离。不要为运行只读检查安装 CUDA，也不要把当前起步包的 NumPy 版本强行用于第三方旧 scGPT 依赖。具体环境锁定在 M0/M4 按真实兼容性完成，见 [environments/README.md](environments/README.md)。

## 已验证代码的范围

`starter/` 包含患者级 manifest 校验、分折、配对效应、同患者四组织对比、组成/状态代数分解、预登记检验与 Holm 校正，以及冻结表示上的透明线性 attention MIL 与 mean 对照。历史 106 项合成/逻辑检查为 33 + 41 + 32，记录见 `starter/evidence/` 及对应结果文件。它不是 scGPT、HA/CTA、MixMIL 的完整复现。

正式 M0 已在隔离的 starter 副本中逐项复核这三组检查。结果与资源记录见 [baseline_checks](artifacts/m0/baseline_checks/)。

## 当前硬件结论

本机 i5-13400 在 CPU 上完成 Pt18 四库的稀疏导入与重读，进程峰值约 662 MiB。M0 探测时项目磁盘可用 480.22 GiB。已枚举到 Intel Arc A750；scGPT 的后端兼容性及资源需求留待 M4 小批实测。

未附大规模原始数据、模型权重及受版权/访问条件约束的其他论文全文。已有审计是可追溯证据，不取代正式实施时的输入核对。论文、第三方代码和数据的原有使用条件仍适用。
