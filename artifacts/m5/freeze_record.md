# V1.1 方法学重构外部验证冻结记录（`CRLM_V1.1_245_F1`）

- **冻结时间戳（UTC）**：`2026-11-21T17:49:04+00:00`
- **冻结版本标识**：`CRLM_V1.1_245_F1`（Git HEAD: `91a589f5ec58736753e5fd71138e0733e3bc8ab1`；初版基线完整封存于 `artifacts/archive/CRLM_V1_245_F1/`）
- **锁定主候选（共 2 项，8 项预登记检验，分属 4 个独立 Holm 校正家族）**：`Mono/Mac/DC::Core_Lipid_SAM_R`（上升）、`T/NK::B07_Th17_TWEAK`（下降）。[冻结候选表](frozen_candidate_registry.csv) · [预登记检验表](registered_tests.csv) · [同机制次要/驻留对照清单](secondary_sensitivity_candidates.csv)
- **方法学与防泄漏保证**：全部候选筛选、20 桶 × 10 重复检测率匹配对照基因集（[matched_controls.csv](matched_controls.csv)）及 `15` 个全开发集迁移分类模型（[frozen_models.csv](frozen_models.csv)）均 100% 在 `GSE245552` 内确定。
- **既往暴露披露**：初版 `CRLM_V1_245_F1` 与纯技术勘误 `CRLM_V1_245_F1_techfix` 对外部队列的历史访问已按秒级 UTC 记录于 [external_access_log.csv](external_access_log.csv)，故 `CRLM_V1.1_245_F1` 严格定位为已披露历史暴露的外部独立队列复现评估。
- **全链路哈希锁**：本记录、验证规则、阶段报告、`scGPT` 预训练权重/词表、全部 `M4` 迁移模型权重及评分/注释/HGNC 源码哈希均已闭环写入 [frozen_inputs_sha256.txt](frozen_inputs_sha256.txt) 与 [output_sha256.txt](output_sha256.txt)。
