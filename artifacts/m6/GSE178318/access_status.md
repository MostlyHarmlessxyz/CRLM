# GSE178318（V1.1 `CRLM_V1.1_245_F1` 次要外部来源状态）

原作者在 GEO 的 Overall design 中明确声明获取或使用表达数据前须联系通讯作者/PI 确认数据使用条件。项目仅保留既有元数据配对与治疗状态审计（共 6 位配对患者，未读取任何细胞表达矩阵），当前尚未获得适用于本研究的书面授权条款（`status = terms_pending`）。

按 `AGENTS.md` 第 6、8 节与 `artifacts/m5/validation_rules.md` 冻结规则，`SECONDARY_DT_178_V1_1` 家族登记的两项次要检验（`SECONDARY_DT_178_V1_1_C1`: `Mono/Mac/DC::Core_Lipid_SAM_R`；`SECONDARY_DT_178_V1_1_C2`: `T/NK::B07_Th17_TWEAK`）原样保留在正式检验清单中：`n_patients` 与 `raw_p` 严格记为缺失（`NA`），仅在家族内 Holm 校正向量中按预登记缺失规则置 `adjustment_only_p = 1.0`（`holm_adjusted_p = 1.0`，`family_size_planned = 2` 不缩减）。该来源在本版本贡献的实际表达验证人数为 0。[原始说明](../../../starter/audit_data/geo/GSE178318_metadata_only.soft) · [登记结果](registered_tests/registered_results.csv)
