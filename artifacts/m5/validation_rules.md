# V1.1 方法学重构版外部验证规则（`CRLM_V1.1_245_F1`）

> [!IMPORTANT]
> **版本属性与历史暴露透明声明**
> `CRLM_V1.1_245_F1` 是在完整归档初版基线 `CRLM_V1_245_F1`（`artifacts/archive/CRLM_V1_245_F1/`）与纯技术勘误 `CRLM_V1_245_F1_techfix` 之后，为修复初版 7 项 CRITICAL 方法学缺陷而在唯一开发集 `GSE245552` 内重做的系统重评估版本。鉴于初版已执行过 `GSE164522` 与 `GSE203608` 评估，本版本在统计与汇报中明确界定为**已披露既往暴露记录的外部独立队列复现与稳健性评估**，不冒称首次盲法确认。

## 1. 冻结主候选与多重比较家族

- `Mono/Mac/DC::Core_Lipid_SAM_R`（预期 `D_T` 与四组织交互 `I` 方向：`上升 (positive)`；焦点纯细类：`Macrophage`）
- `T/NK::B07_Th17_TWEAK`（预期 `D_T` 与四组织交互 `I` 方向：`下降 (negative)`；焦点纯细类：`CD4 T cell`）

- **四个独立校正家族（每族双侧精确符号检验 + Holm，`alpha = 0.05`）**：
  1. `PRIMARY_DT_164_V1_1`：`GSE164522` 患者级配对肿瘤差 `D_T = LM - PT`；
  2. `BACKGROUND_I_164_V1_1`：`GSE164522` 同患者四组织合格子集 `Q` 上的器官背景扣除交互 `I = D_T_Q - D_N_Q`；
  3. `SECONDARY_DT_178_V1_1`：`GSE178318` 次要方向证据（使用条款待核实前保留 `raw_p = NA`，仅 Holm 向量置 `1.0`）；
  4. `SECONDARY_DT_203608_V1_1`：`GSE203608` 次要方向与精度证据。

## 2. 统一基因轴、冻结背景扣除与细胞清洗规则

1. **HGNC 规范化共同核基因轴**：所有外部队列统一经 `src/crlm/hgnc_alias_map.py` 规范化符号（`GSE203608` 结合 `ENSG` 消歧重名符号），并在与 `GSE245552` 共享的核编码基因轴（`GSE164522` 采用 `245_164_shared_nuclear` `18,815` 核基因；`GSE203608` 采用 `245_164_203608_shared_nuclear` `18,733` 核基因，显式排除全部 `353` 个 MT/核糖体偏倚基因）上评分。冻结正向成员测量覆盖须 `>= 80%` 且 `>= 5` 个基因。
2. **主评分公式（零外部重拟合）**：主指标固定为 `debiased_ucell_scores`（全测量轴 Mann–Whitney UCell `max_rank=1500` 减去本冻结包锁定的 `GSE245552` 20 桶 × 10 重复检测率匹配非成员对照集 `matched_controls.csv`），**严禁在外部队列重算检测率分桶或重抽背景基因**。同步输出共同核基因重标定 `log1p_cp10k_debiased`、`ucell_raw` 与初版 `legacy_positive` 敏感性对照。
3. **跨队列一致细胞清洗与同规则重注释对照**：
   - `GSE164522`：主分支默认剔除原作者已标注的全部跨谱系 `CD3E+` 双细胞群（髓系 `hC56/hC41/hC47/hC36` 共 3,662 个细胞 + B/Plasma `hC75_B-CD3E` 3,174 个细胞），另附保留双细胞敏感性与 **`annotation_m1.classify` 统一分层标志物重注释对照分支**；
   - `GSE203608`：严格过滤体内 `Sample_pass_qc == TRUE`、`Cell_type != "NA"` 且 `MitoFraction <= 0.20`（剔除 12,136 个死细胞/未注释细胞），同时提供洁净作者标签主分支与 `annotation_m1.classify` 统一分层标志物重注释分支，并附剔除患者 49 敏感性。
4. **聚合与分解**：每个合格标本×谱系门槛 `>= 50` 个细胞（另附 `20/100` 及纯细类 `>= 20` 对照），按“文库细胞均值 → 同患者同组织病灶等权 → 患者配对差”汇总；组成比例统一报告 `total_resolved_immune` 与 `core_three_lineage_immune` 双分母，细类状态执行共同支持对称分解。
