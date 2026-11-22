# GSE178318（V1.1 `CRLM_V1.1_245_F1` 严格前瞻性零前序表达暴露次要外部验证状态）

当前状态：`status = authorized_prospective_zero_exposure_validated`（已获正式使用授权并完成冻结后严格前瞻性零前序表达暴露验证）。

## 1. 前瞻性零前序表达暴露密码学证据链

- **冻结先行（零前序表达/效应暴露）**：在 `CRLM_V1.1_245_F1` 预登记冻结（`2026-11-21T17:49:04+00:00`）及本地代码与模型哈希归档提交（Git commit `5b01467`）完成前，项目仅审阅公开论文（Che et al., *Cell Discovery* 2021）与 GEO 样本配对/治疗状态/条形码元数据，**从未下载或读取 `GSE178318` 任何单细胞表达矩阵字节，亦未计算任何基因程序或分类器效应**（`artifacts/m5/frozen_inputs_sha256.txt` SHA256 = `5ce11557b7bc565184ad7418befc08fec0fe411ab11586de7a92a8d1c7314f14`；`artifacts/m5/output_sha256.txt` SHA256 = `0b52cc300ebaece035572287ee71144d31a5a789e5a51e0435468f19ec1682fe`）。
- **授权后解锁与原封不动执行**：用户在 `5b01467` 提交后正式获得 `GSE178318` 数据使用授权并指令解锁执行。`scripts/fetch_gse178318.py` 在首次下载首字节前首先核对 `artifacts/m5/{frozen_inputs_sha256,output_sha256}.txt` 与提交 `5b01467` 字节级一致（加固后的现行脚本与全仓校验同步逐项核验全部 `72 + 19` 个底层冻结文件的 SHA256 零偏移），随后于 `2026-11-22T03:05:15Z–03:21:35Z` 下载并校验 `barcodes.tsv.gz`、`genes.tsv.gz`、`matrix.mtx.gz`；全过程对 `artifacts/m5/`（含 `frozen_inputs_sha256.txt` 锁定的全部 72 个路径）、`configs/project.json`、`src/crlm/{annotation_m1,program_scores,hgnc_alias_map,model_bags}.py` 及 `models/m4/full/*` 保持 **0 字节修改**。

## 2. 样本配对、PBMC 物理隔离与正式检验完成度

- 全套 15 个文库（`140,281` 个条形码，来自 6 位同步切除 IV 期 MSS 结直肠癌肝转移患者：`COL07, COL12, COL16` 未新辅助治疗 vs `COL15, COL17, COL18` 术前化疗）已完成逐条形码核验（见 [15 文库审计表](coverage/all_15_gsm_barcode_audit.csv)）。
- 其中 3 个外周血文库（`GSM5387664` `COL12_PBMC`、`GSM5387671` `COL17_PBMC`、`GSM5387674` `COL18_PBMC`，共 `26,950` 个细胞）已按协议 **100% 物理隔离排除**，仅将 12 个配对实体肿瘤文库（`6 PT + 6 LM = 113,331` 个输入细胞，零改动 M2 QC 后保留 `105,813` 个细胞）物化至 `data/processed/m6/GSE178318/*.h5ad`。
- `SECONDARY_DT_178_V1_1` 家族两项预登记检验均在完整 **N = 6** 位配对患者中实测完成（两项均符合预登记方向 `matches_registered_direction`），使全研究预登记检验表实现 **8 / 8 项 100% 实测闭环（0 项缺失）**。[源输入清单](source_inputs.csv) · [外部访问日志](external_access_log.csv) · [登记检验结果](registered_tests/registered_results.csv) · [阶段报告](stage_report.md)
