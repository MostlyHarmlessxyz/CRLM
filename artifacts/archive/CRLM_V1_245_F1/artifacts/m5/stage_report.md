# M5 · 候选及验证冻结

已从 GSE245552 的开发结果锁定髓系 Mac_SPP1 上升和 B/Plasma UPR 下降。冻结于 2026-11-21 12:44:41 UTC，对应开发版本 `ff1f1a3`。当时已登记外部队列的身份、测量范围和既往发表亚群暴露。新增外部程序效应在此之后读取。[冻结记录](freeze_record.md) · [既往暴露](external_access_log.csv)

`scripts/freeze_validation.py` 已在 M4 Python 环境实际运行。`starter/adjust_registered_tests.py` 的计划验证、两个唯一候选、八个唯一检验、四个各含两项的家族均已核对。两个主候选的正向测量成员、四个模型程序的成员、共同排序轴、单细胞评分代码和 15 个全开发分类头可逐项定位。运行 `sha256sum -c artifacts/m5/frozen_inputs_sha256.txt` 已通过。[主要候选](frozen_candidate_registry.csv) · [全部检验](registered_tests.csv) · [规则](validation_rules.md) · [模型](frozen_models.csv)
