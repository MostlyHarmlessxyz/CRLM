# 阶段执行与复算约定（`workflow/`）

项目按 `M0 → M1 → M2 → M3 → M4 → M5（密码学冻结）→ M6 → M7 → M8` 顺序执行，阶段门禁定义见 [IMPLEMENTATION_PLAN.md](../docs/IMPLEMENTATION_PLAN.md)，完整复算命令序列与 26 项脚本入口映射见 [reproduce.md](../artifacts/m8/reproduce.md) 与 [reproducibility_manifest.csv](../artifacts/m8/reproducibility_manifest.csv)。

## 核心规则

1. **阶段可追溯性**：每个阶段在 `artifacts/m0/` … `artifacts/m8/` 下输出正式 `stage_report.md`、统计结果表与 `output_sha256.txt` 清单。
2. **单向防泄漏**：`M5` 冻结后不得根据 `M6` 外部队列结果反向修改 `V1.1` 候选定义、特征映射或分类头权重（由 `artifacts/m5/frozen_inputs_sha256.txt` 强制核验）。
3. **存储分层**：原始矩阵、细胞表示缓存与模型权重分别存放于 `data/`、`embeddings/`、`models/`，不混入轻量阶段产物目录。
