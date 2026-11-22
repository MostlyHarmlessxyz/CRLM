# GSE335552 多尺度无损空间邻近分析规则

1. **材料与独立样本单位**：公开 Xenium 数据集 `GSE335552` 含 4 张原位切片（`M242T`、`M131T`、`M397T` 肝转移灶及 `M131N` 匹配邻近肝），对应 **3 位独立患者**（无原发结直肠组织，故仅支持肿瘤内细胞来源与局部空间关系，不构成 PT–LM 配对复现）。
2. **`CRLM_V1.1_245_F1` 对齐面板审计**：使用 `src/crlm/hgnc_alias_map.py` 结合探针符号与 `ENSG` 标识规范化基因轴；对齐 `artifacts/m5/gene_sets/` 与 `artifacts/m3/program_definitions_draft/program_genes.csv`，在 `B19_S4_Mac_SPP1_clean` 中同时剥离跨谱系 T/NK 颗粒毒力基因（`GZMB`、`GNLY`、`LAT`、`IL2RG`）与内皮污染基因 `FLT1`（全集 `159` 基因 / 共享核轴 `157` 基因），并同步核查 `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`（全集 `113` / 共享核轴 `111` 基因）、`Core_Lipid_SAM_R`（`16` 基因）、`B07_Th17_TWEAK`（`14` 基因）、`Resident_Kupffer`（`14` 基因）及 `Core_Chemokine_IFN_C`（`14` 基因）探针覆盖率。
3. **接触边界转录本溢出消歧（双分支对照）与 `~EPCAM` 边界说明**：
   - 原始候选：`CD8 T` 要求 `(CD3D|CD3E) & (CD8A|CD8B) & ~EPCAM`；髓系要求 `CD68/CD163/C1QA/CD14/LYZ` 至少检出 2 项、含前 3 项之一且 `~EPCAM`。
   - **主分支（`dominance_disambiguated`）**：当重叠细胞满足 `CD3D+CD3E+CD8A+CD8B > CD68+CD163+C1QA+CD14+LYZ` 时保留为接触态 `CD8 T`；满足髓系转录本 `>= 2 ×` T 细胞转录本时保留为接触态髓系；仅剔除近等量模糊重叠细胞，两个集合严格互斥。
   - **敏感性分支（`strict_zero_spillover`）**：剔除全部跨谱系重叠细胞作透明对照。
   - **`~EPCAM` 巢缘溢出注意事项**：紧贴高表达 `EPCAM+` 肿瘤上皮巢（`<15 µm`）的真实巨噬细胞若捕获 `>= 1` 条溢出的 `EPCAM` 转录本，会被 `~EPCAM` 规则排除，从而减少 `<15 µm` 近巢层可观测髓系分母并引入巢缘截断选择，解读近巢 vs 间质层时须与网格平移检验一并审视。
4. **多尺度网格、半格平移（`+25 µm` / `0.5 × tile`）边界核查与精确置换检验**：同时评估 `400 / 200 / 100 / 50 µm` 四档方格及 `SPP1 >= 1 / >= 2` 门槛，拆分 `all_myeloid`、`near_epcam_lt15um`（`<15 µm`）与 `stromal_epcam_ge15um`（`>=15 µm`），并附加坐标原点平移 `+25 µm`（及各尺度半格 `0.5 × tile`）敏感性核查（`grid_shift_sensitivity.csv`）。所有置换采用 `BASE_SEED = 20261121` 经 SHA-256 独立派生的随机种子执行 `1000` 次格内标签置换，报告双侧精确 Phipson–Smyth `p = (b + 1) / 1001`（实测耗时 `26.21s`，峰值 RSS `610.6 MiB`）。
