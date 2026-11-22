# M3 · GSE245552 患者级程序与组成

**Mono/Mac/DC 的 Mac_SPP1 状态标记在全部 11 位合格患者的 LM 中高于 PT。** 平均配对差为 +0.01974，患者 bootstrap 95% 区间为 +0.01067 至 +0.03094。B/Plasma 的未折叠蛋白应答（UPR）也通过预定开发筛选，平均差为 −0.01219，7/8 位患者下降；其平均差区间包含 0，且对细胞测序深度分布敏感。两项均为有待独立复现的开发结果。Mac_SPP1 采用已发表的巨噬细胞状态基因块。

## 患者效应

处理 39 个 GSM 的 190,031 个 M2 保留细胞。每 GSM 每谱系至少有 50 个细胞才能计算主状态读数，组成分母包括所有保留细胞。相同患者、相同组织的合格病灶均值等权。下表为患者内 LM−PT 的平均差，分数按四种组织共同提供的 29,569 个基因进行单细胞内排序。区间按患者有放回重抽样 10,000 次得到。

| 谱系、文献读数 | 配对患者 | 平均差及 95% 区间 | 同向患者 | 去掉任一患者时的均差范围 |
|---|---:|---|---:|---|
| Mono/Mac/DC · Mac_SPP1 | 11 | +0.01974（+0.01067 至 +0.03094） | 11/11 | +0.01530 至 +0.02150 |
| B/Plasma · UPR | 8 | −0.01219（−0.02355 至 +0.00339） | 7/8 | −0.01909 至 −0.00998 |
| Mono/Mac/DC · Mac_CXCL10 | 11 | −0.00571（−0.01960 至 +0.00921） | 7/11 | −0.01021 至 −0.00239 |
| T/NK · IFN-γ 应答 | 12 | +0.00138（−0.00727 至 +0.01027） | 6/12 | −0.00094 至 +0.00384 |

输入材料提供各程序的原始成员，包括 164、113、94 和 200 个基因；245 共同测量轴上分别有 160、107、90 和 196 个。Mac_CXCL10 与 T/NK IFN-γ 应答均未达到预定的患者同向比例。原程序主得分、辅助通路交集、两组非肿瘤背景及全部覆盖状况可从源表重建。[开发规则与原始引用](program_definitions_draft/development_rules.md) · [基因与覆盖](program_definitions_draft/program_genes.csv) · [逐库主读数](sample_program_scores.csv) · [逐患者差值](patient_effects/patient_contrasts.csv) · [全部效应与区间](patient_effects/effects_summary.csv)

**四组织。** 在 Mac_SPP1 的唯一完整患者 Pt18 中，D_T_Q 为 +0.01828，D_N_Q（NL−NC）为 −0.01167，I 为 +0.02995。B/Plasma 虽有一个邻近组织配对，Pt18 的 LM 中该谱系少于 50 个细胞，故其四组织共同子集为空。两组织和四组织的记录保留各自患者分母。

## 细胞组成、深度与病灶

13 位有肿瘤配对的患者均可估计捕获细胞组成。LM 相对 PT 的 B/Plasma 比例平均下降 15.40 个百分点，12/13 下降；Mono/Mac/DC 平均上升 3.54 个百分点，10/13 上升。这些值适用于经过 CD45 磁选重混的样本。程序分解在每个父谱系内使用 M2 的身份子组。共同子组每侧至少 20 个细胞时，Mac_SPP1 的类型内状态分量中位数为 +0.0153，组成分量中位数为 +0.0002；对应 PT、LM 共同支持质量中位数分别为 0.855、0.978。UPR 两种分量的中位数分别为 −0.0044、−0.0054。逐人的分量可加，中位数各自汇总。[组成差](composition_results/effects_summary.csv) · [共同支持分解](composition_results/patient_subtype_decomposition.csv)

Mac_SPP1 去掉其他被评分程序集的共享成员后均差为 +0.01765，100 细胞门槛下为 +0.01873（8 对）。10 组检测率匹配随机基因的平均差从 −0.00130 至 +0.00360。UPR 去共享后为 −0.01174，而 10 组随机集合的平均差也都为负，从 −0.00582 至 −0.00030。UPR 的患者配对差与该谱系每细胞可见计数中位数之差相关，Spearman ρ=0.833（8 对）。[基因与覆盖敏感性](patient_effects/sensitivity_summary.csv) · [配对深度相关](paired_depth_sensitivity.csv)

针对该深度相关，在发现集对四个完整程序统一补做探索性标准化。按单细胞可见总计数的 1,000、3,000、10,000 切分四档，将两侧各有至少 10 个细胞的分档按两组织比例的较小值重标化，保持病灶等权。Mac_SPP1 的标准化均差为 +0.01833，10/11 位仍上升；UPR 缩至 −0.00400，5/8 位仍下降，两侧深度分布可匹配质量的中位数为 0.733。深度同时反映真实 RNA 含量和测序因素，这种标准化只能描述共同深度范围的敏感性。UPR 适合作为捕获 B/Plasma 群体的弱候选。[分档后的逐患者效应](patient_depth_standardized.csv) · [开发短名单](candidate_evidence_draft.md)

## 伪批量基因差异

按病灶级 ≥50 同谱系细胞选样，在每个谱系内只纳入 PT、LM 两侧都合格的患者。对作者已筛选的整数 counts 按患者、组织求和，拟合 PyDESeq2 0.5.3 负二项模型，设计为 `~ donor_uid + tissue`，对照为 LM 相对 PT。各谱系独立做 BH。多病灶病例另用每个患者每种组织的首个合格标本重算。六次拟合均通过样本顺序和设计秩验证。

| 谱系 | 患者对 | 合并病灶的受检基因 | BH q < 0.05 | 首个病灶受检基因 | 首个病灶 BH q < 0.05 |
|---|---:|---:|---:|---:|---:|
| T/NK | 12 | 10,407 | 841 | 10,209 | 648 |
| Mono/Mac/DC | 11 | 9,055 | 1,405 | 8,893 | 1,313 |
| B/Plasma | 8 | 9,292 | 47 | 9,030 | 71 |

髓系 SPP1 的 log₂ 倍数变化为 +4.56，q=8.15×10⁻⁸；首病灶分析为 +4.55，q=4.19×10⁻⁷。B/Plasma 的 XBP1 和 HSPA5 在合并病灶分析中 q 分别为 0.995、0.722。髓系、B/Plasma 的参数离散度趋势拟合由软件切换到均值趋势。逐基因收敛状态在结果中保留，B/Plasma 合并病灶的 IGKV2-29 LFC 未收敛（q=0.247）。独立计数分支提供该捕获混合内的基因证据，核心程序估计仍采用病灶等权秩分数。[六次拟合概览](paired_counts/fit_summary.csv) · [参考基因](paired_counts/reference_readouts.csv) · [全部基因、设计及文库](paired_counts/)

## 重建与交接

从项目根目录运行。M2 的 H5AD、共同基因轴及伪批量索引是输入。单细胞环境沿用 M2，负二项模型使用单独锁定的环境。程序与组成分析耗时 38.90 秒，峰值 RSS 501 MiB。配对负二项六次拟合耗时 64.85 秒，父进程峰值约 288 MiB。另用 M0 样本表直接重算 160 个患者×程序条目的差值，核对了 390 条样本分数的 M2 细胞数、84 行共同支持分解恒等式及六次 NB 拟合的设计秩和结果基因轴。[秩评分资源、种子与输入 SHA](resource_profile.json) · [NB 资源](paired_counts/profile.json) · [环境锁](paired_counts_packages.txt)

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 uv run --no-project --python .venv-m1/bin/python python -Bu scripts/analyze_discovery_programs.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 uv run --no-project --python .venv-m1/bin/python python -Bu scripts/analyze_paired_depth.py
de_env="$HOME/.local/share/venvs/crlm-m3-de"
uv venv --python 3.12.13 "$de_env"
uv pip install --python "$de_env/bin/python" -r artifacts/m3/paired_counts_packages.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 uv run --no-project --python "$de_env/bin/python" python -Bu scripts/analyze_paired_counts.py
```

40 组匹配随机读数仅用发现集估计检测率，属于技术对照。M4 内部比较可使用四个外部文献预定义完整集合，筛选通过的两个程序作为 M5 候选草案。外部效应与确认结果归 M6。[秩评分代码](../../src/crlm/program_scores.py) · [数据与输出索引 SHA](output_sha256.txt)
