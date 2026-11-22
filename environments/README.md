# 仓内运行时环境与依赖锁清单（`environments/`）

本目录保存项目四个隔离 Python 3.12 虚拟环境的精确 `uv` 锁文件，确保跨阶段数值计算与模型推理可复现。

## 1. 环境分工

| 环境目录 | 锁文件 | 核心依赖 | 负责阶段 |
|---|---|---|---|
| `.venv-m1` | `venv-m1.lock.txt` | `numpy`, `pandas`, `scipy`, `anndata`, `h5py` | **M0–M1**：`starter/` 106 项自检、`GSE245552` 原始稀疏导入与审计 |
| `.venv-m2` | `venv-m2.lock.txt` | `.venv-m1` + `matplotlib`, `pillow` | **M2**：全开发集自适应质控、两阶段免疫注释与组成分解 |
| `.venv-m4` | `venv-m4.lock.txt` | `scikit-learn`, `statsmodels`, `pydeseq2`, `pyarrow`, `openpyxl` | **M3–M8 主线**：HGNC 对齐、去偏 UCell、配对伪批量、PCA/程序分类器、M5 冻结、外部程序打分、Xenium 空间分析与回归测试 |
| `.venv-scgpt` | `venv-scgpt.lock.txt` | `torch==2.15.0+xpu`（Intel Arc A750 XPU）, `scikit-learn`, `statsmodels` | **M4 & M6 深度表征**：冻结 `scGPT` 全细胞编码、`HA/CTA/MixMIL` 对照训练与外部队列模型推理 |

> 四个环境统一使用 `Python 3.12.13` 及相同版本的 `numpy==2.5.5`、`pandas==3.1.1`、`scipy==1.18.1`、`anndata==0.13.4`。`scripts/check_project.py` 仅依赖标准库。

## 2. 从锁文件重建环境

```bash
uv venv --python 3.12 .venv-m1 && uv pip install --python .venv-m1/bin/python -r environments/venv-m1.lock.txt
uv venv --python 3.12 .venv-m2 && uv pip install --python .venv-m2/bin/python -r environments/venv-m2.lock.txt
uv venv --python 3.12 .venv-m4 && uv pip install --python .venv-m4/bin/python -r environments/venv-m4.lock.txt
uv venv --python 3.12 .venv-scgpt && uv pip install --python .venv-scgpt/bin/python \
  --extra-index-url https://pytorch-extension.intel.com/release-whl/stable/xpu/us/ \
  -r environments/venv-scgpt.lock.txt
```

## 3. 快速导入自检

```bash
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -B -c "import statsmodels, patsy, pydeseq2, openpyxl, sklearn, anndata, pyarrow"
PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -B -c "import torch, statsmodels, sklearn; assert torch.xpu.is_available()"
```
