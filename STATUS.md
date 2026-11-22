# CRLM 项目状态台账（`CRLM_V1.1_AUDIT_REMEDIATED`）

## 版本与验证属性

- **当前主版本（`CRLM_V1.1_245_F1`）**：已完成 `M0–M8` 全链路分析与 8/8 项预注册外部检验。其中 `GSE164522` 与 `GSE203608` 因初版曾读取表达矩阵，界定为**已披露前序暴露的外部稳健性重评估**；`GSE178318`（6 对同步切除 IV 期 MSS 患者）在 M5 密码学锁定后首次解锁原始 UMI 矩阵并零改动运行，构成**零表达暴露外部验证**（跨队列未经基因型去重，供者非重叠状态保持 `unproven`）。
- **历史对照归档**：初版不可变快照封存于 [artifacts/archive/CRLM_V1_245_F1/](artifacts/archive/CRLM_V1_245_F1/)；纯客观技术勘误轨存于 `artifacts/techfix_v1_f1/`。

## 阶段完成总表（M0–M8）

| 阶段 | 状态 | 核心交付与结论摘要 | 阶段报告 |
|---|---|---|---|
| **M0** 环境与输入审计 | completed | 核实 `GSE245552`（16 人 / 13 对 PT–LM / 39 库），通过 106 项合成基线自检，完成官方 HGNC 快照、`scGPT` 语料重叠与外部暴露审计。 | [M0 报告](artifacts/m0/stage_report.md) |
| **M1** 最小真实导入 | completed | 完成 `Pt18` 四组织（22,179 细胞）稀疏导入、源坐标核验与两阶段无冲突免疫层级注释。 | [M1 报告](artifacts/m1/stage_report.md) |
| **M2** 全开发集处理 | completed | 39 库自适应质控保留 173,450 个细胞（已解析免疫谱系内 `unresolved = 0%`）；`>=50` 细胞门槛下 `T/NK`、`Mono/Mac/DC`、`B/Plasma` 分别含 12、11、8 对配对患者。 | [M2 报告](artifacts/m2/stage_report.md) |
| **M3** 患者程序与配对计数 | completed | 20 桶检出率匹配去偏 UCell + 供者共识 cNMF + 配对伪批量（Student-$t$ 与 $2^N$ 符号翻转）：锁定脂质/SAM 巨噬细胞（11/11 升）与黏膜 Th17–TWEAK（10/12 降）两项主候选；证伪淘汰 `B/Plasma::UPR`（纯浆细胞伪批量 0 DEG）。 | [M3 报告](artifacts/m3/stage_report.md) |
| **M4** 全袋分类与 MIL 基准 | completed | A750 XPU 编码 70,628 个肿瘤免疫细胞；完成折内线性 PCA vs 零样本 `scGPT`、均值池化 vs 注意力 MIL、`HA/CTA`、`MixMIL`、200 次置换零分布、`256×10` 子袋方差审计及 15 个全集迁移头拟合。 | [M4 报告](artifacts/m4/stage_report.md) |
| **M5** 原子预注册冻结 | completed | 在解锁 `GSE178318` 表达矩阵前锁定 2 个主候选、8 项预注册检验（4 个 Holm 家族）、15 个迁移头与 72 项上游输入 SHA-256（[冻结记录](artifacts/m5/freeze_record.md)）。 | [M5 报告](artifacts/m5/stage_report.md) |
| **M6** 外部复现与盲测验证 | completed | `GSE164522` 宽谱系 `T/NK` 四联体观测交互显著通过家族校正（$I = -0.049$，Holm $P = 0.043$），纯 CD4+ T 细胞在 `164`（9/10 降）与盲测 `178`（5/6 降）一致下调；脂质/SAM 巨噬细胞在盲测 `178` 达 5/6 升（宽谱系未校正肿瘤对比 0/2 过 Holm）；双外部队列 15/15 迁移头 100% HVG 覆盖零拒判，线性 PCA 整体匹配或优于零样本 `scGPT`。 | [M6 报告](artifacts/m6/stage_report.md) |
| **M7** 空间审计与验证设计 | completed（湿实验 0 项） | `GSE335552`（4 张 10x Xenium 切片）救回 3,584 个接触态 CD8+ T 细胞，`400 → 50 µm` 多尺度网格与 `<15 µm` 肿瘤巢分层量化中尺度构型效应；制定四组织盲法 `SPP1+CXCL9` 原位杂交 + 自体髓系条件培养基 CD8+ T 细胞 `Transwell` 趋化方案（含 `NT_Ctrl / SPP1_KD / SPP1_Rescue_Endo / rhSPP1_Ext` 组别；已完成实际湿实验数为 0）。 | [M7 报告](artifacts/m7/stage_report.md) |
| **M8** 终稿集成与复算归档 | completed | 交付三版本对照终稿、出版级主图与逐元素源数据表、局限偏差记录及 26 项全核实复算清单（[最终报告](artifacts/m8/final_report.md) · [复算指南](artifacts/m8/reproduce.md)）。 | [M8 报告](artifacts/m8/stage_report.md) |

## 防篡改治理规则

`configs/project.json`、`artifacts/m5/`、`models/m4/full/` 与 `src/crlm/` 严格保持 `GSE178318` 解锁前的 M5 锁定哈希不变（由 `sha256sum -c artifacts/m5/frozen_inputs_sha256.txt` 强制核验）；`GSE225857` 因缺失公开患者配对映射维持排除状态。
