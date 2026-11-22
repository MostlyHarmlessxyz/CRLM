# 全链路确定性复算指南（`CRLM_V1.1_AUDIT_REMEDIATED` & 三版本归档）

所有命令均从仓库根目录 `/home/mostly_harmless/Projects/CRLM` 执行，统一建议前缀 `PYTHONDONTWRITEBYTECODE=1`（并限制 `OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4` 控制多核 BLAS 内存峰值）。各阶段脚本、核验产物表、SHA256 清单、隔离 `.venv-*` 解释器及依赖锁文件逐行登记于 [reproducibility_manifest.csv](reproducibility_manifest.csv)。

## 一、环境映射与非破坏性门禁自检

| 隔离虚拟环境 | 依赖锁清单 | 负责阶段与用途 |
|---|---|---|
| `.venv-m1/bin/python` | `artifacts/m1/python-packages.txt` | M1 `Pt18` 四组织稀疏矩阵导入与源坐标核对 |
| `.venv-m2/bin/python` | `artifacts/m2/python-packages.txt` | M0 结构契约检查、M2 全开发集 39 库导入与两阶段层级注释 |
| `.venv-m4/bin/python` | `artifacts/m4/python-packages.txt`（附 `artifacts/m3/paired_counts_packages.txt`） | M3 去偏程序/配对伪批量、M4 PCA/scGPT 发现集基准与迁移头拟合、M5 冻结、M6 外部程序验证与分类推理（`predict`）/汇总、M7 空间审计、Track 1 程序勘误与模型评估、M8 主图渲染 |
| `.venv-scgpt/bin/python` | `artifacts/m4/resource_profiles/xpu-packages.txt`（附 `artifacts/m4/direct_prior/extra_packages.txt`） | Intel Arc A750 XPU（**`torch 2.15.0+xpu`**）72 细胞精度预检、M4/M6/Track 1 全袋 `scGPT` CLS 编码，以及 M4/Track 1 先验 MIL（`HA/CTA/MixMIL`）评估 |

快速结构、硬件与全仓回归自检（不修改任何受控仓库文件，仅在系统临时目录运行 106 项隔离合成检查）：

```bash
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/check_project.py
PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python scripts/preflight.py --torch-probe
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/test_audit_remediations.py
```

---

## 二、`CRLM_V1.1_AUDIT_REMEDIATED`（主交付轨）端到端复算命令

```bash
# 1. M1–M2：开发集 GSE245552 稀疏导入、官方 HGNC 对齐、库自适应 QC 与两阶段零冲突层级注释（173,450 保留细胞）
PYTHONDONTWRITEBYTECODE=1 .venv-m1/bin/python -Bu scripts/import_gse245552.py --manifest data/manifests/m0/sample_manifest_validated.csv --donor GSE245552:Pt18
PYTHONDONTWRITEBYTECODE=1 .venv-m2/bin/python -Bu scripts/process_gse245552.py

# 2. M3：19 个完整程序 20 桶检测率匹配去偏 UCell 筛选、供者共识 cNMF、三层次配对负二项（Student-t + 符号翻转）与深度敏感性
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/analyze_discovery_programs.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/analyze_paired_counts.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/analyze_paired_depth.py

# 3. M4：A750 XPU 72 细胞预检 + 全袋 70,628 肿瘤免疫细胞编码、去混杂消融、200 次置换零分布、256x10 子袋方差、内层调参先验 MIL 与 15 个全集迁移头
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/encode_scgpt.py --profile --device xpu
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/encode_scgpt.py --device xpu --batch 16
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/benchmark_discovery.py --representation PCA
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/benchmark_discovery.py --representation scGPT
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/benchmark_prior_mil.py --representation PCA32 --device xpu
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/benchmark_prior_mil.py --representation scGPT32 --device xpu
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/benchmark_prior_mil.py --representation scGPT512 --device xpu
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/summarize_benchmarks.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/fit_transfer_models.py

# 4. M5：直接校验已锁定的 CRLM_V1.1_245_F1 72 项全链路输入哈希（注：若加 --force 重跑 freeze_validation.py 会写入新的 UTC 冻结时间戳）
sha256sum -c --quiet artifacts/m5/frozen_inputs_sha256.txt
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/freeze_validation.py --self-check

# 5. M6：GSE164522（四组织 + 洁净/同规则双分支）、GSE203608 洁净体内验证、全袋外部编码与 15/15 零拒判分类推理及三版本对照汇总
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse164522.py identity
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse164522.py convert --freeze-commit SHA_LOCKED
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse164522.py score --freeze-commit SHA_LOCKED
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse203608.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/predict_gse164522.py development-check
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/predict_gse164522.py --freeze-commit WORKTREE --device xpu --batch 16 --overwrite encode
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/predict_gse164522.py --freeze-commit WORKTREE predict
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/assemble_validation.py

# 6. M7–M8：GSE335552 多尺度/肿瘤巢空间邻近审计与四幅出版级主图（PNG + 矢量 PDF + 7 份源表）重建
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/spatial_support.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/build_paper_figures.py
```

---

## 三、原版归档核验与 Track 1（`CRLM_V1_245_F1_techfix`）纯技术勘误复算

```bash
# 核验不可变原版基线归档（CRLM_V1_245_F1）
(cd artifacts/archive/CRLM_V1_245_F1 && sha256sum -c --quiet archive_sha256.txt)
sha256sum -c --quiet models/archive/CRLM_V1_245_F1/ARCHIVE_SHA256.txt
sha256sum -c --quiet embeddings/archive/CRLM_V1_245_F1/ARCHIVE_SHA256.txt

# 复算 Track 1 纯客观 Bug 修复分步消融（锁定初版候选、方向、评分公式与冻结模型权重）
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/run_techfix_v1_f1_programs.py
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/run_techfix_v1_f1_models.py mixmil-bugfix
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/run_techfix_v1_f1_models.py --device xpu --batch 4 encode-external
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/run_techfix_v1_f1_models.py evaluate-external
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/run_techfix_v1_f1_models.py finalize
```

---

## 四、冷启动原始公开数据拉取参考（`GSE245552`、`GSE164522`、`GSE203608`、`GSE335552`）

- `GSE245552`：`PYTHONDONTWRITEBYTECODE=1 .venv-m2/bin/python scripts/fetch_gse245552.py`
- `GSE164522`：`PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/validate_gse164522.py fetch`
- `GSE203608`：`PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/fetch_gse203608.py`
- `GSE335552`（仅拉取本研究使用的 8 个细胞特征矩阵与质心表，跳过超大原始图像归档）：

```bash
mkdir -p data/raw/GSE335552
awk '/^!Sample_supplementary_file = / && /(_cell_feature_matrix.h5|_cells.parquet.gz)$/ {sub("^ftp:", "https:", $3); print $3}' starter/audit_data/geo/GSE335552_metadata_only.soft |
while IFS= read -r url; do
  curl -fLsS --retry 3 -o "data/raw/GSE335552/${url##*/}" "$url" || exit 1
done
```
