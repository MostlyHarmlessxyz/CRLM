# 本地与云端算力指南

**针对配置：Intel Core i5-13400、64 GB RAM、约 400 GB 可用 SSD、Intel Arc A750。核验日期：2026-11-21。** 操作系统、驱动、现有 Python/R 环境和目标机实际可用内存尚未现场检查；以下为工程建议，运行时间和峰值占用不是本机实测结果。

## 1. 当前决定：本地启动，暂不购买云套餐

这台机器适合完成项目的大部分常规计算：文献与数据审计、按样本稀疏导入、可计算的 QC、注释/程序评分、患者级配对分析、counts 分支的伪批量统计、PCA/程序/组成基线，以及冻结表示上的小型 MIL。前提是采用分样本、分谱系和有内存预算的处理，避免一次致密化所有队列。

**M0–M1 从本地 CPU 开始；不需要等 GPU 环境准备好，也不需要先购买 AutoDL/Colab。** M2–M3 可继续本地推进。到 M4 scGPT 编码阶段，根据实际 profile 决定继续 CPU、试验 A750 XPU，或短租 NVIDIA GPU。先做一次可靠的冻结表示编码，再把表示带回本地比较下游模型，是当前合理分工。

“本地完成”在这里指计算研究工作。新增患者组织、空间测量和湿实验需要相应实验条件；算力不能替代这些证据。Codex 选择 GPT-6.2 Astra 与给 scGPT 租 GPU 是两件事：Agent 帮助编写/运行项目，科研模型运行在你选择的计算环境。这里不建议或提供把 Astra 权重部署到 A750 的方案。

## 2. 各工作适合在哪里做

| 工作 | 默认位置 | 原因与约束 |
|---|---|---|
| 文献、manifest、样本身份与来源审计 | 本地 CPU | 几乎不涉及 GPU，正确身份比加速优先 |
| M1 最小真实稀疏导入、轴核对与 QC 能力 | 本地 CPU | 按样本读取；基因缺失不能靠大内存修复 |
| M2 全开发集、分谱系处理与注释 | 本地 CPU 优先 | 64 GB 可用来规划分块；全队列联合图/多副本的占用须先测 |
| UCell/固定程序、患者表、bootstrap/精确检验 | 本地 CPU | 按细胞块和基因范围计算，最后以患者作为独立单位 |
| 合法 counts 的伪批量统计 | 本地 CPU/R | 聚合后规模较小；164 的连续表达仍不能套 counts 模型 |
| 程序/组成、PCA mean、透明 attention MIL | 本地 CPU | 起步包本来就是 CPU 实现；限制嵌入同时驻留及调参规模 |
| scGPT 冻结编码 | 先 profile，必要时短租 NVIDIA | 序列长度、batch、依赖和中间激活决定资源，不能只看参数文件大小 |
| 完整 HA/CTA、MixMIL 等直接前作复现 | 按实际模型 profile | 简单基线先完成；不把所有方法同时排入多进程扫描 |
| scGPT 全模型微调、多 GPU 或从头预训练 | 当前主线不安排 | 本项目患者数与研究目标不支持把这作为起步必要条件 |
| 实测空间面板检查与部分空间分析 | 本地按样本 | 先拿小的细胞×基因表；原始大图像只在所需阶段取用 |

这张表是根据方案工作量给出的实施判断，不是某台机器跑完全部数据的性能认证。可按 `artifacts/m1/resource_profile.json` 与 M4 profile 调整。

## 3. A750：可以研究 XPU，不能假定原版 scGPT 直接可用

Intel 官方标注 A750 为 **8 GB GDDR6**。PyTorch 当前 Intel GPU 指南把 Arc A-Series 列入支持/验证范围，使用 `xpu` 后端，并要求匹配的系统和驱动。官方也提供 `torch.xpu.is_available()` 检查。因此“A750 完全不能做深度学习”不准确；但框架级支持没有验证本项目锁定的 scGPT 及其全部依赖。[C1,C2]

本报告审计的 scGPT commit 为 `cebd6fae655b9c585a4807daa3ac31bb764f06b4`。该版本 `cell_emb.py` 默认 `device="cuda"`，编码路径使用 `torch.cuda.amp.autocast`，默认 fast-transformer 指向 flash；还会在构造细胞数据集前对输入执行 `toarray()`。README 说明 flash-attention 可选并存在 CPU 路线，但这不等于原实现已支持 XPU。[C3,C4]

因此 A750 路线需要单独确认设备、dtype/autocast、注意力实现、DataLoader、必要算子及依赖。只把字符串改为 `"xpu"` 不足以证明可用。新版 PyTorch XPU 与旧第三方依赖的版本约束也可能冲突。应将试验放在独立环境，保留原 CPU 路径。

建议顺序：

1. M0 先记录 OS、Python/R、驱动和已安装包；无 PyTorch 也能完成轻量探测及其余 M0 工作。
2. 如已有 PyTorch，执行 `python scripts/preflight.py --torch-probe`，只确认是否识别后端。这个命令不运行 scGPT，也不安装软件。
3. 只有进入 M4 且准备试用 A750 时，依当时的 [PyTorch Intel GPU 指南](https://docs.pytorch.org/docs/stable/notes/get_start_xpu.html) 建立隔离环境。不要机械安装 NVIDIA CUDA toolkit，也不要把某个旧 IPEX/oneAPI 教程当作当前唯一方案。
4. 用一小块已核实细胞完成前处理、编码和结果回写；检查形状、有限数值、细胞顺序及缓存 provenance。固定同一基因 token 子集，比较 CPU 与 XPU 路径的数值和表示相似性，容差依据 dtype 设定并记录，不能要求混合精度逐位相同。
5. 观察设备使用、实际吞吐和峰值显存，排除大量 CPU fallback 后再决定扩大。若移植变成独立工程问题，停止扩大范围，改用短租 CUDA 环境完成已定义的编码任务。

**不承诺 A750 比 CPU 更快，也不承诺 8 GB 对任意序列长度/batch 足够。** A750 试验不是科学主线的先决条件。

## 4. 64 GB 内存：主要风险是数据表示与复制

历史 245 审计记录为 190,190 个输入细胞（尚未独立 QC），肿瘤共同基因轴为 29,569。以整个细胞集合取共同基因轴构造单份致密矩阵作示例：

\[
190190\times29569\times4\;\mathrm{bytes}\approx22.50\;\mathrm{GB}\approx20.95\;\mathrm{GiB}.
\]

float64 则约 45.0 GB / 41.9 GiB。以上只是一份数值矩阵，未计 AnnData 元数据、raw/layers、归一化副本、算法工作数组和其他队列。**64 GB 并不能保证把所有队列转为 dense 后仍然稳定。** 输入数量来源见完整报告第 3 章和 `starter/audit_data/`；算式可由 `scripts/preflight.py` 直接重算。[C6]

推荐按样本和谱系保留稀疏矩阵，在需要时分块读取；查看第三方函数内部是否转 dense。scGPT 要在调用其内层编码前就切块，不能只减 GPU batch 却把全输入留给外层 `toarray()`。同时固定逐细胞 token 抽样种子，使改变分块顺序不悄悄改变表示；记录截断和缺基因比例。[C3]

冻结表示的量级小很多。例如，若实际所选 checkpoint 的表示维度为 512，190,190 个细胞的 float32 向量本体约 **389.5 MB**；100 万个细胞约 **2.048 GB**。这是示例维度和裸数组大小，非所选模型已经确定、也不含索引和多版本副本。这解释了为何“云上编码一次，回本地做下游对照”有实际价值。

并发初值建议为 **一个重计算作业、四个 worker**，之后由实际 CPU/I/O/内存调整。这是保守起步值，不是对 i5-13400 的最优性能结论；避免 Agent 同时运行多份全量导入或嵌入编码。

## 5. 400 GB SSD：可用于当前路线，按阶段取数

400 GB 是开始时空闲量，不是可以一次全部用完的下载额度。规划时保留约 **80–100 GB** 给系统、临时解压、环境和意外峰值；这是建议余量，不是已经测得的项目固定消耗。

先下载 M1 样本所需的 MTX、features、barcodes 和元数据；M2 再扩到开发集。不要先下全部 FASTQ、全部图谱或空间原始图像。245 的 39 个压缩 MTX 历史审计合计约 1.056 GB，但这不含所有外部输入、解压文本、派生对象和环境，不能据此把全项目预算写为 1 GB。[C6]

每次新增大对象，记录“原下载 + 解压暂存 + 最终稀疏对象/缓存 + 保留旧版本”在同一时刻的占用。对确定可重建且不需保留的中间物采用明确清理规则；原始输入与最终冻结成果仍应保留可靠副本。模型权重按实际 checkpoint 获取一份并校验，不为比较后端重复下载多套相同权重。

已有空间审计也显示：GSE335552 先用约 35 MB 的四份细胞×基因 H5 就能检查面板，不能因为发布页有约 26 GB 包而默认全部下载。是否需要原始影像，取决于 M7 的具体问题。[C6]

## 6. 什么时候值得租 GPU

**现在不租。** M4 前先得到：确定的输入细胞与基因范围、固定 checkpoint/词表、明确的前处理和 token 抽样、实际小批量 profile。以下任一条件成立，才有理由转云端：

- CPU 全量预计时长超过你为编码设定的可接受等待时间；估计基于代表性小批次，不凭型号猜小时数。
- A750 需要的后端移植和依赖排障超出你愿意投入的工程范围。
- 实际显存/内存不足，且在保持预定算法和信息公平性的条件下，减 batch/分块仍不能解决。
- 后续确有已定义、资源测量支持的模型比较需要 GPU，而不是为了“以防万一”先买月租。

供首轮试租的工程起点是 **单张 NVIDIA GPU、16–24 GB 显存**；预算允许时 24 GB 留有更多调整空间。这个范围是待 profile 的选择建议，**不是所有 scGPT 配置可运行的官方最低/充分规格**。需要同时检查实例 CPU RAM 与磁盘，不能只比较 GPU 名称。当前没有证据要求先买 A100/H100、多卡或包月。

## 7. AutoDL 与 Colab 的选择

| 选项 | 本项目适用场景 | 需要接受的取舍 |
|---|---|---|
| 本地 | 全部准备、常规统计、调试、报告、部分基线；结果和数据主副本 | scGPT 全量 CPU 编码可能较慢，A750 兼容未测 |
| AutoDL 按量实例 | 已确定输入和环境后的正式批量编码；需要持续保留环境、明确选 GPU 的作业 | 自行管理环境和结果备份，按实例开机时间计费，型号/价格以实时页面为准 |
| Colab 免费/个人付费 | 小型 Notebook、官方示例、一次短 profile 或教学试验 | GPU 类型、可用额度、运行时长和中断条件动态，不适合作为本项目唯一的长时无人值守主环境 |

AutoDL 官方计费文档说明：按量实例按**开机到关机**计费，不是只按 GPU 忙碌时间；关机停止实例计算收费，但付费扩容数据盘可能继续计费。因而最好先在本地准备输入和运行脚本，上云后完成编码、校验并下载结果，再关机；释放实例前先确认成果已有可读副本。[C5]

Colab 官方 FAQ 明确资源不是无限或保证供应，GPU 类型、可用性和运行时限制会变化；付费个人方案也存在动态资源条件。需要按批保存结果和断点，不能用单次长 Notebook 作为全部研究的唯一记录。[C7]

**若到 M4 确需购买，我优先建议短租 AutoDL 一台合适的 NVIDIA 实例；已有稳定 Colab 使用经验时，也可用 Colab 做同样的小批验证。** 这是针对本项目需要可控环境和批处理的工程判断，不是平台性能实测排名。当前不用同时买两家，也不预报固定人民币总价。

费用估计应记录：

\[
\text{估计总费用}=\text{实例开机小时}\times\text{当时小时单价}+\text{存储及其他明确费用}.
\]

开机时间包含安装、下载、调试、编码和回传期间仍开机的部分。完整吞吐预测还要分别计入预处理/读写、模型加载、预热与稳定编码；小批外推只给条件化范围。不要把 GPU kernel 时间当账单时间。

## 8. 最小 profile 记录与上云回传

M1 记录稀疏导入的真实样本/细胞/基因数、矩阵 dtype/nnz、耗时、峰值 RSS、最终文件和临时目录体积。M4 再记录模型 commit、权重/词表 SHA、backend、dtype、max length、实际 batch、token 抽样策略、设备、峰值显存、吞吐、OOM/警告和输入/输出哈希。

profile 可先用少量细胞检查正确性，再用更大但受限的一批检查吞吐是否稳定。具体规模由谱系、内存和初次结果决定，不把固定“1,000 个细胞”当作一定足够的科学测试。CPU/XPU/CUDA 的比较使用相同预定输入与 token，不允许因为某后端更方便而无记录地改变基因范围或模型。

如使用云端，带回的最小成果为：按 sample/细胞 ID 索引的 embeddings、同版本 manifest、配置、权重和词表哈希、代码 commit、环境锁、执行日志、profile、输出校验和。先本地读回并核对细胞顺序、维度、有限数值和重复/缺失，再关闭或释放资源。云端失败或后端更换不自动改变外部验证角色。

## 9. 当前即可做的操作

1. 把完整工程解压到本地常用项目目录；保留 400 GB 空间的阶段预算。
2. 用 Codex Astra 的 Plan 模式读取本项目，要求明确 M0–M1 输入、输出与验收。
3. 进入执行阶段，在本机运行 `scripts/check_project.py` 和 `scripts/preflight.py`。结果只代表这台运行机器。
4. 复核原 starter 的目标环境，再做最小真实稀疏导入。GPU 尚未配置不会阻断这一步。
5. 到真正需要 scGPT 编码时再比较 CPU profile、A750 试验成本与短租 NVIDIA 的费用。届时再买算力。

## 官方来源与本项目证据

- **[C1]** [Intel Arc A750 官方规格](https://www.intel.com/content/www/us/en/products/sku/227954/intel-arc-a750-graphics/specifications.html)：8 GB GDDR6。
- **[C2]** [PyTorch：Getting Started on Intel GPU](https://docs.pytorch.org/docs/stable/notes/get_start_xpu.html)：当前支持范围、系统/驱动要求、XPU 安装及验证。实际安装当天重新核对。
- **[C3]** [scGPT 固定 commit 的 cell_emb.py](https://github.com/bowang-lab/scGPT/blob/cebd6fae655b9c585a4807daa3ac31bb764f06b4/scgpt/tasks/cell_emb.py)：默认 CUDA、CUDA autocast、fast backend 和致密化路径；原文件哈希见 `starter/source_audit/`。
- **[C4]** [scGPT 官方 README，固定 commit](https://github.com/bowang-lab/scGPT/blob/cebd6fae655b9c585a4807daa3ac31bb764f06b4/README.md)：可选 flash-attention 与 CPU/GPU 加载说明；不把旧 CUDA 安装示例直接当作当前 Intel 安装指南。
- **[C5]** [AutoDL 官方计费](https://api.autodl.com/docs/price/)；[快速开始](https://www.autodl.com/docs/quick_start/)：实例、时长和存储计费以购买时页面为准。
- **[C6]** [完整研究报告第 3 章](CRLM_Research_Plan.md)、`starter/audit_data/` 的原始元数据/矩阵审计：输入细胞、gene 集、压缩文件与空间面板体积；本指南资源算式为透明估算，不是训练结果。
- **[C7]** [Google Colab 官方 FAQ](https://research.google.com/colaboratory/faq.html)：动态资源、GPU 类型、运行时及保存限制。
