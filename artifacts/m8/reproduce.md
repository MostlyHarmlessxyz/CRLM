# 全链路确定性复算指南（`CRLM_V1.1_AUDIT_REMEDIATED` & 三版本归档）

所有命令均从仓库根目录 `/home/mostly_harmless/Projects/CRLM` 执行，统一建议前缀 `PYTHONDONTWRITEBYTECODE=1`（并限制 `OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4` 控制多核 BLAS 内存峰值）。各阶段直接生产脚本、代表性直接产物表、SHA256 清单、收归仓内的可执行 `.venv-*` 解释器及依赖锁文件逐行登记于 [reproducibility_manifest.csv](reproducibility_manifest.csv)（共 26 条全核实记录；注：早期 M0 预检与 `starter` 清单核验曾使用系统 `/usr/bin/python` 及 `~/.local/share/venvs/crlm-m0-starter`，见 [M0 报告](../m0/stage_report.md)，现已全部收归仓内 `.venv-m1`/`.venv-m2` 统一复算）。

## 一、环境映射与非破坏性门禁自检

| 隔离虚拟环境 | 依赖锁清单 | 负责阶段与用途 |
|---|---|---|
| `.venv-m1/bin/python` | `artifacts/m1/python-packages.txt` | M1 `Pt18` 四组织稀疏矩阵导入与源坐标核对 |
| `.venv-m2/bin/python` | `artifacts/m2/python-packages.txt` | M0 供者清单校验、M2 全开发集 39 库导入与两阶段层级注释 |
| `.venv-m4/bin/python` | `artifacts/m4/python-packages.txt`（附 `artifacts/m3/paired_counts_packages.txt`） | M3 去偏程序/配对伪批量、M4 发现集 PCA/scGPT 评估与迁移头拟合、M5 冻结、M6 外部程序验证与 `164` `predict`/汇总、M7 空间审计、Track 1 程序与外部模型评估、M8 主图渲染 |
| `.venv-scgpt/bin/python` | `artifacts/m4/resource_profiles/xpu-packages.txt`（附 `artifacts/m4/direct_prior/extra_packages.txt`） | Intel Arc A750 XPU（**`torch 2.15.0+xpu`**）72 细胞精度预检、M4/`164`/`178`/Track 1 全袋与子袋 `scGPT` CLS 编码、`178` 前瞻分类推理，以及 M4/Track 1 先验 MIL（`HA/CTA/MixMIL`）评估 |

快速结构、硬件与全仓回归自检（不修改任何受控仓库文件；其中 `starter` 106 项合成检查复制至系统临时目录 `$TMPDIR` 隔离运行以保持 `starter/` 只读）：

```bash
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/check_project.py
PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python scripts/preflight.py --torch-probe
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/test_audit_remediations.py
```

---

## 二、`CRLM_V1.1_AUDIT_REMEDIATED` 与 Track 1 端到端因果依赖复算命令

> [!NOTE]
> **M5 密码学锁保护约定**：`artifacts/m3/`、`artifacts/m4/{full_data_models.csv,output_sha256.txt,stage_report.md}`、`artifacts/m5/`、`models/m4/full/`、`configs/project.json`、`src/crlm/` 与 `starter/` 共 72 项路径已在 `5b01467` 提交中由 `artifacts/m5/frozen_inputs_sha256.txt` 锁定。若需重跑步骤 1–3（尤其原地重跑 M4 会更新耗时 JSON 的运行时间戳从而偏离已锁定 M5 哈希），请在隔离工作副本或独立输出目录执行；常规外部复核直接从**步骤 4** 起跑即可。

```bash
# 1. M1–M2：开发集 GSE245552 稀疏导入、官方 HGNC 对齐、库自适应 QC 与两阶段零冲突层级注释（173,450 保留细胞）
PYTHONDONTWRITEBYTECODE=1 .venv-m1/bin/python -Bu scripts/import_gse245552.py --manifest data/manifests/m0/sample_manifest_validated.csv --donor GSE245552:Pt18
PYTHONDONTWRITEBYTECODE=1 .venv-m2/bin/python -Bu scripts/process_gse245552.py

# 2. M3：19 个完整程序 20 桶检测率匹配去偏 UCell 筛选、供者共识 cNMF、三层次配对负二项（Student-t + 2^N 符号翻转）与深度敏感性
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

# 4. M5 与原版归档校验：直接核验已锁定的 CRLM_V1.1_245_F1 72 项全链路输入哈希及不可变 V1 归档（不带 --force 重跑以免覆盖冻结时间戳）
sha256sum -c --quiet artifacts/m5/frozen_inputs_sha256.txt
sha256sum -c --quiet artifacts/m5/output_sha256.txt
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/freeze_validation.py --self-check
(cd artifacts/archive/CRLM_V1_245_F1 && sha256sum -c --quiet archive_sha256.txt)
sha256sum -c --quiet models/archive/CRLM_V1_245_F1/ARCHIVE_SHA256.txt
sha256sum -c --quiet embeddings/archive/CRLM_V1_245_F1/ARCHIVE_SHA256.txt

# 5a. M6 单队列程序验证与双外部队列全袋 XPU 编码 / 15/15 分类迁移（含 GSE178318 真实三版本子袋 XPU 编码）
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse164522.py identity
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse164522.py convert --freeze-commit SHA_LOCKED
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse164522.py score --freeze-commit 5b01467
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse203608.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/validate_gse178318.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/predict_gse164522.py development-check
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/predict_gse164522.py --freeze-commit WORKTREE --device xpu --batch 16 --overwrite encode
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/predict_gse164522.py --freeze-commit WORKTREE predict
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/predict_gse178318.py --freeze-commit 5b01467

# 5b. Track 1（CRLM_V1_245_F1_techfix）纯客观技术勘误分步消融（在 assemble_validation 与出图前运行以供三版本汇总读取）
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/run_techfix_v1_f1_programs.py
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/run_techfix_v1_f1_models.py mixmil-bugfix
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 .venv-scgpt/bin/python -Bu scripts/run_techfix_v1_f1_models.py --device xpu --batch 4 encode-external
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/run_techfix_v1_f1_models.py evaluate-external
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/run_techfix_v1_f1_models.py finalize

# 5c. M6 家族 Holm 校正、8/8 全观测检验与三版本并列对照总表组装
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/assemble_validation.py

# 6. M7–M8：GSE335552 多尺度/半格平移/肿瘤巢空间邻近审计与四幅确定性出版级主图（PNG + 矢量 PDF + 9 份逐元素源表）重建
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/spatial_support.py
PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python -Bu scripts/build_paper_figures.py
```

---

## 三、冷启动原始公开数据拉取参考（`GSE245552`、`GSE164522`、`GSE178318`、`GSE203608`、`GSE335552`）

- `GSE245552`：`PYTHONDONTWRITEBYTECODE=1 .venv-m2/bin/python scripts/fetch_gse245552.py`
- `GSE164522`：`PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/validate_gse164522.py fetch`
- `GSE178318`（确认具备原文献数据使用许可后执行；命中本地 SHA256 时幂等保留原 `downloaded_utc`）：`PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/fetch_gse178318.py`
- `GSE203608`：`PYTHONDONTWRITEBYTECODE=1 .venv-m4/bin/python scripts/fetch_gse203608.py`
- `GSE335552`（仅拉取本研究使用的 8 个细胞特征矩阵与质心表，跳过超大原始图像归档）：

```bash
mkdir -p data/raw/GSE335552
awk '/^!Sample_supplementary_file = / && /(_cell_feature_matrix.h5|_cells.parquet.gz)$/ {sub("^ftp:", "https:", $3); print $3}' starter/audit_data/geo/GSE335552_metadata_only.soft |
while IFS= read -r url; do
  curl -fLsS --retry 3 -o "data/raw/GSE335552/${url##*/}" "$url" || exit 1
done
```
