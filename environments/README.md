# 环境策略

本交接包的 `scripts/` 仅用 Python 标准库，不需先安装分析依赖。正式执行时再按 M0 锁定操作系统和 Python/R 版本。

| 环境 | 用途 | 当前状态 |
|---|---|---|
| `crlm-starter` | NumPy/Pandas/SciPy 的原起步代码 | Python 3.12.14 + `starter/requirements.txt` 是历史已验证组合；目标主机待复核 |
| `crlm-scrna` | 稀疏导入、Scanpy/AnnData、注释、评分 | 尚未实现/锁定，M0–M2 按真实输入建立 |
| `crlm-rstats` | 需要整数 counts 的伪批量统计等 R 工作 | 尚未实现/锁定；不对 GSE164522 处理后连续值套 counts 模型 |
| `crlm-scgpt` | 固定权重编码、CPU/CUDA/XPU 后端试验 | M4 单独建立；先审计依赖及硬件，不能强行沿用起步包 NumPy 版本 |

不提供未经安装验证的巨型 `requirements.txt` 或自动 GPU 安装脚本。创建后记录精确包版本、平台、依赖锁、所用来源/commit、权重/词表哈希和安装日志。环境与项目目录分开；数据、环境缓存和下载压缩包均计入磁盘预算。

先用现有 Python 运行轻量 preflight。原起步代码若要复跑，在独立虚拟环境中使用其 `requirements.txt`。本次不自动创建虚拟环境或安装依赖。
