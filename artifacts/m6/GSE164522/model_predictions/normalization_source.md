# GSE164522 外部肿瘤分类输入与共同核基因轴重标定（Track 2 · `CRLM_V1.1`）

1. **原始归一化来源与逐标本尺度核验**：
   - Liu 等（*Cancer Cell* 2022, STAR Methods e3, “Single-cell RNA-Seq data processing”）说明 10x Genomics 5′ 液滴数据经质控后按每细胞文库总量归一到 `10,000`（`CP10k`），再施加 `scanpy.pp.log1p()`。公开肿瘤表达表的 `74,888` 个细胞键与原 metadata 100% 唯一匹配（[论文](https://doi.org/10.1016/j.ccell.2022.02.013) · [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164522)）。
   - 在作者完整特征轴上逐标本核验非空细胞 `sum(expm1(X))`：中位数严格落入 `9,900–10,100`，且 `>= 95%` 细胞位于 `9,500–10,500`（[PT 尺度](../coverage/PT_scale_check.csv) · [LM 尺度](../coverage/MT_scale_check.csv)）。

2. **跨平台共同核编码基因轴对齐与 `CP10k` 核子轴重标定**：
   - 由于 `GSE245552` 肿瘤矩阵缺失 353 个正常组织特有特征（含全部 13 个线粒体基因及数百个高丰度核糖体蛋白基因），直接使用含核糖体/MT 基因分母的原始 `CP10k` 会系统性压低核基因幅度（约 20%–35% 分母偏移）。
   - `V1.1` 统一经 `src/crlm/hgnc_alias_map.py` 将 `GSE164522` 特征无碰撞映射至共同核编码基因轴 `245_164_shared_nuclear`（**`18,815 / 18,815 = 100.0%` 全量命中**），在每个细胞内经 `expm1` 恢复线性比例后按核基因总和重新归一至 `10,000` 并取 `log1p`，使外部输入尺度与 `GSE245552` 训练输入严格同构。

3. **跨谱系双细胞剔除、全袋保留与零拒判特征覆盖**：
   - 在建袋前剔除原作者已注释的全部跨谱系 `CD3E+` 髓系双细胞（`hC56/hC41/hC47/hC36`）与 `hC75_B-CD3E` B–T 双细胞，并保留每个合格标本×谱系（`>= 50` 个细胞）的**全部洁净细胞**（`max_cells=None`，消除 256 细胞抽样噪声）。
   - 因 `M4` 全集 `PCA` 头严格在 `245_164_shared_nuclear` 上训练，三谱系全部训练 HVG 在 `GSE164522` 均实现 **`512 / 512 = 100.0%` 精确覆盖（0 缺列、0 中心替代、0 模型拒判）**（[HVG 测量记录](pca_measurement.csv)）。
   - `program_elastic` 直接调用 `models/m4/full/<lineage>/` 锁定的 `program_defs.json` 与 `matched_controls.csv` 计算深度去偏 `debiased_ucell_scores`（零外部背景重拟合）；`scGPT` 在同一重标定核基因轴上经 `18,717` 个单射词表 token 执行 51-bin 分位数分箱、`max_length=1200` 逐细胞确定性截断与 L2 归一化。
