# M0 通过

本机复核得到 39 个公开 GSM、16 位患者、13 对 PT–LM 和 1 个四组织患者。Pt18 可用同一组四份标本进入 M1。真实来源、配对计数和队列级协议见 [manifest 核查](manifest_review.md) 与 [运行机验证后的样本表](../../data/manifests/m0/sample_manifest_validated.csv)。

按 starter 默认的 5 折和固定种子 20261121，已生成 [QC 前患者分折](../../data/manifests/m0/donor_folds_pre_qc.json)。各折留出 4、3、3、3、3 位患者，39 个 GSM 各留出一次，同人全部病灶与组织归入同折。此分配只使用身份表；建模时再按 M2 各谱系的实际覆盖核算每折人数。

| 本机合成检查 | 通过项 | 耗时（秒） | 最大 RSS（KiB） | 退出码 |
|---|---:|---:|---:|---:|
| 统计与分割 | 33 | 7.692 | 151364 | 0 |
| 线性 attention 与均值对照 | 41 | 6.211 | 140024 | 0 |
| 登记检验与 Holm | 32 | 10.074 | 148608 | 0 |

三组共 106 项，本机结果逐项与历史检查对应，原始三份 JSON 的 SHA 保持一致。每次运行的 UTC 时间与进程资源记录见 [run_metrics.csv](baseline_checks/run_metrics.csv)，原始日志和结果 JSON 在 [隔离副本](baseline_checks/starter/)。合成检查沿用各自固定种子。M0 初始主机探测为 2026-11-21 11:43:21 UTC；最后一项登记自检于 11:44:04 UTC 启动，用时 10.074 秒。

## 重建

从项目根目录运行。Bash 示例中的检查目录应使用单独的本次输出路径。

```bash
root="$PWD"
envdir="$HOME/.local/share/venvs/crlm-m0-starter"
uv run --no-project --python /usr/bin/python python -B scripts/preflight.py --output artifacts/m0/environment_probe.json
uv venv --python 3.12.13 "$envdir"
uv pip install --python "$envdir/bin/python" -r starter/requirements.txt
uv run --no-project --python "$envdir/bin/python" python -B starter/validate_manifest.py --manifest starter/audit_data/GSE245552_sample_manifest.csv --outdir data/manifests/m0
uv run --no-project --python "$envdir/bin/python" python -B starter/make_splits.py --donor-manifest data/manifests/m0/donor_manifest_validated.csv --n-folds 5 --seed 20261121 --output data/manifests/m0/donor_folds_pre_qc.json
checks="$root/artifacts/m0/baseline_checks/starter"
mkdir -p "$checks"
cp starter/*.py "$checks/"
cp -R starter/synthetic_demo/. "$checks/synthetic_demo/"
cd "$checks"
export PYTHONDONTWRITEBYTECODE=1
for script in selfcheck.py mil_selfcheck.py registered_tests_selfcheck.py; do
  uv run --no-project --python "$envdir/bin/python" python -B "$script" > "$script.log" 2>&1 || exit 1
done
```

本次在同一 uv 环境内用 `wait4` 记录三项的退出码和子进程峰值 RSS。逐 GSM 依据 GEO 标题、来源组织、患者和 URL 与两个样本表及逐库 JSON 连接；对照 S1 的患者临床字段、S2 的患者×组织数及 S3 的作者组与末列 gene。

## 证据与范围

输入、源码 SHA 与所用历史结果见 [输入清单](input_inventory.csv)，阶段产物 SHA 见 [输出清单](output_sha256.txt)。本地共有 29 份论文 PDF，32 条本地归档记录。85 条主参考的阅读深度沿用原表。初始环境和已辨识设备见 [运行环境](environment_report.md)。

外部队列沿用 [角色决定](data_role_decisions.md) 及 [暴露日志](external_access_log.csv)。临床逐人时序、178 使用条件、225 患者映射及原作者 PhenoAligner 实现各自留待对应分析。M0 的真实数据工作止于已归档来源和身份核查。生物学效应、真实训练与候选冻结仍待相应阶段。
