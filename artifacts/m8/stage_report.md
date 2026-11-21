# M8 · 本地 V1 交付

正式报告附四幅主图，每幅提供 PNG、矢量 PDF 与源数据。患者程序、同一四组织子集、模型覆盖和空间读数分别呈现；两项主要校正检验及两条受限来源记录全部保留。[最终报告](final_report.md) · [图](figures/) · [源数据](source_data/) · [限制与未测项](limitations_and_deviations.md)

阶段索引覆盖 M0–M8 的源码、数据表、环境及运行记录。冻结版本为 `72ee86b`，外部联合分析为 `72a0f5d`。原始矩阵、逐细胞对象、嵌入和冻结模型在项目约定的数据目录内；轻量交付包含原件定位、版本、完整结果和重建入口。[复算命令](reproduce.md) · [阶段索引](reproducibility_manifest.csv)

已从正式阶段表重建全部主图，重读图源核对患者数、同一四组织差分、八项登记方向和族内 Holm 结果。164 的独立逐细胞复算与程序实现的最大差为 5.6×10⁻¹⁷，分类预测核对了 315 个方法×标本×谱系位置和缺失原因。四幅图均进行了版面核对。阶段交付文件和图源可用下列命令核验。[外部单元验证](../m6/GSE164522/acceptance.json) · [外部分类验证](../m6/GSE164522/model_predictions/stage_report.md)

```bash
sha256sum -c artifacts/m8/output_sha256.txt
```

本地可用来源的 V1 分析与实验准备已交付。GSE178318 的 PI 条款、GSE225857 的供者身份、未获得的受控资源及实际人源干预均有明确边界。新来源或实验数据取得后可以单独扩展研究。[来源与纳入状态](../m6/cross_cohort_evidence_cards.md) · [实验材料](../m7/material_and_readout_requirements.csv)
