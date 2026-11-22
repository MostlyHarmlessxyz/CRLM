# 仓内运行时环境与依赖锁清单（`environments/`）

本目录记录 CRLM 项目全部阶段（`M0–M8`、`Track 1 techfix` 与 `Track 2 V1.1`）唯一受控的四个仓内 Python 虚拟环境精确锁文件，彻底消除仓外隐式环境（如 `~/.local/share/venvs/`）与临时 `PYTHONPATH` 注入（如 `~/.cache/crlm-m4-priors/pylibs`）。

## 1. 主机与解释器基线

- **操作系统**：Linux x86_64
- **硬件配置**：Intel Core i5-13400（CPU 主线）、64 GB RAM、Intel Arc A750 Discrete GPU（`xpu` 后端专用于冻结 `scGPT` 编码与 PyTorch MIL 对照）
- **Python 解释器**：`Python 3.12.13`（四个仓内 `.venv-*` 完全统一）
- **包管理器**：`uv`（所有锁文件均由 `uv pip freeze --python <venv>/bin/python` 导出）
- **共享数值底座**：四个环境统一锁定 `numpy==2.5.5`、`pandas==3.1.1`、`scipy==1.18.1`、`anndata==0.13.4`、`h5py==3.16.0`、`zarr==3.4.2`，确保跨阶段稀疏矩阵、`.h5ad` 与统计结果二进制读写零版本漂移。

## 2. 仓内环境与锁文件对照表

| 仓内环境路径 | 锁文件 | 包数量 | 核心增量依赖 | 唯一负责的项目阶段与脚本 |
|---|---|---:|---|---|
| `.venv-m1` | `environments/venv-m1.lock.txt` | 27 | `numpy`, `pandas`, `scipy`, `anndata`, `h5py` | **M0 基线自检 & M1 稀疏导入**：`starter/*selfcheck.py`、`scripts/stage_gse245552.py`、`scripts/audit_gse245552.py` |
| `.venv-m2` | `environments/venv-m2.lock.txt` | 34 | 在 `.venv-m1` 基础上增加 `matplotlib==3.11.2`, `pillow==13.0.0` | **M1–M2 质量控制、分层注释与组成分解**：`scripts/process_gse245552.py`、M2 诊断绘图 |
| `.venv-m4` | `environments/venv-m4.lock.txt` | 61 | `scikit-learn==1.9.1`, `statsmodels==0.15.0`, `patsy==1.0.3`, `pydeseq2==0.6.0`, `openpyxl==3.1.5`, `pyarrow==26.0.1`, `matplotlib==3.11.2` | **M0 HGNC 审计、M3 程序发现与伪批量 `PyDESeq2` 统计、M4 透明线性 MIL/PCA 基准、M5 候选冻结、M6 外部验证打分（`GSE164522`/`GSE203608`）、M7 Xenium Parquet 空间分析、M8 图表与回归测试**：`scripts/analyze_discovery_programs.py`、`scripts/analyze_paired_counts.py`、`scripts/benchmark_discovery.py`、`scripts/freeze_validation.py`、`scripts/validate_gse164522.py`、`scripts/validate_gse203608.py`、`scripts/assemble_validation.py`、`scripts/spatial_support.py`、`scripts/build_paper_figures.py`、`scripts/test_audit_remediations.py` |
| `.venv-scgpt` | `environments/venv-scgpt.lock.txt` | 70 | `torch==2.15.0+xpu`（Intel SYCL/oneMKL XPU 运行时）, `scikit-learn==1.9.1`, `statsmodels==0.15.0`, `patsy==1.0.3` | **M4 & M6 冻结 `scGPT` 细胞表示编码（Intel Arc A750 XPU）及 PyTorch 直接前作 MIL 对照（`HA`/`CTA`/`MixMIL`，直接读取仓内 `references/code/{hier-mil,mixmil}`）**：`scripts/encode_scgpt.py`、`scripts/benchmark_prior_mil.py`、`scripts/predict_gse164522.py` |

> **标准库脚本说明**：`scripts/check_project.py` 与 `scripts/run_m0_preflight.py` 仅依赖 Python 标准库，可使用系统 `python3` 或任一 `.venv-*/bin/python` 运行。

## 3. 从锁文件重建环境的标准命令

```bash
uv venv --python 3.12 .venv-m1 && uv pip install --python .venv-m1/bin/python -r environments/venv-m1.lock.txt
uv venv --python 3.12 .venv-m2 && uv pip install --python .venv-m2/bin/python -r environments/venv-m2.lock.txt
uv venv --python 3.12 .venv-m4 && uv pip install --python .venv-m4/bin/python -r environments/venv-m4.lock.txt
# .venv-scgpt 需配置 Intel PyTorch XPU wheel 索引以解析 torch==2.15.0+xpu：
uv venv --python 3.12 .venv-scgpt && uv pip install --python .venv-scgpt/bin/python --extra-index-url https://pytorch-extension.intel.com/release-whl/stable/xpu/us/ -r environments/venv-scgpt.lock.txt
```

## 4. 环境完整性快速自检

```bash
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -B -c "import statsmodels, patsy, pydeseq2, openpyxl, sklearn, anndata, pyarrow"
PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -B -c "import torch, statsmodels, sklearn; assert torch.xpu.is_available()"
```
