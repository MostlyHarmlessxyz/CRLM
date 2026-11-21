# CRLM 项目状态 · 2026-11-21

## 已完成的准备工作

- 完整研究设计、数据与文献审计，以及对应 Markdown 报告。
- 原起步包 106 项合成/逻辑检查的历史记录。此记录不能代替目标机器的验证。
- Codex 交接骨架、阶段合同、Plan/执行提示词、论文与算力准备清单。
- 本次交接包的实际可用性核验仅写入 `artifacts/handoff/`。

## 正式研究阶段

| 阶段 | 状态 | 已有证据或待交付结果 |
|---|---|---|
| M0 环境与输入 | completed | 39 个 GSM、16 位患者、13 个肿瘤配对；本机 106 项合成检查通过。[M0 报告](artifacts/m0/stage_report.md) |
| M1 最小真实导入 | completed | Pt18 四组织，1 位患者、22,179 个细胞；4 个稀疏对象重读、原始计数抽查与初级谱系覆盖。[M1 报告](artifacts/m1/stage_report.md) |
| M2 全开发集 | completed | 39 库共 190,190 个细胞，统一 QC 保留 190,031 个。每病灶至少 50 个相应细胞时，T/NK、Mono/Mac/DC、B/Plasma 的 PT–LM 患者对分别为 12、11、8。[M2 验收](artifacts/m2/stage_report.md) |
| M3 患者程序 | completed | 髓系 Mac_SPP1 在 11/11 对中升高；B/Plasma UPR 在 7/8 对中降低，深度标准化后为 5/8。六个配对负二项拟合及组成/状态分解已归档。[M3 报告](artifacts/m3/stage_report.md) |
| M4 模型比较 | completed | 390 次拟合，2,279 个折外袋预测；在相同患者与袋上完成 PCA/scGPT、均值/attention 和三个直接前作；15 个迁移头已拟合。[M4 报告](artifacts/m4/stage_report.md) |
| M5 冻结 | pending | 程序定义、方向、范围、检验家族与验证规则锁定 |
| M6 外部复现 | pending | 独立来源效应、同一共同子集背景与失败记录 |
| M7 空间/功能 | pending | 实测空间支持、可实施实验及其证据边界 |
| M8 交付 | pending | 可复现研究结果、图表、方法与全部局限 |

GSE245552 的状态估计基于作者交付且可测的基因与固定谱系标记。组成分母是各库捕获并通过 QC 的细胞。真实候选仍为开发草案。

## 更新规则

执行者只有在实际输出与日志满足 `docs/IMPLEMENTATION_PLAN.md` 的门槛后，才能同步更新本文件和 `configs/project.json` 的阶段状态。不得把脚本存在、计划写成或合成示例通过等同于真实分析完成。保留失败、不可测以及外部暴露记录。
