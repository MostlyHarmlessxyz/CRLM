# M4 · 层级 MIL 与 MixMIL 直接对照

相同输入和患者折下，三种直接前作均未显示跨谱系一致的 AUROC 增益。PCA32 髓系的 MixMIL 为 0.994，均值头为 0.988，差 +0.006，按患者重抽样的 95% 区间为 0.000 至 +0.031；T/NK 两者均为 1.000。原生 scGPT512 的层级头和降至 32 维后的三个直接头在各谱系均低于对应均值头。评估目标为已取得组织的来源，生物复现证据另按 M5–M6 的患者程序检验取得。

| 输入及模型 | T/NK | Mono/Mac/DC | B/Plasma |
|---|---:|---:|---:|
| PCA32 均值 | 1.000 | 0.988 | 0.867 |
| PCA32 CTA | 0.759 | 0.970 | 0.658 |
| PCA32 HA | 0.933 | 0.970 | 0.617 |
| PCA32 MixMIL | 1.000 | 0.994 | 0.833 |
| scGPT32 均值 | 0.918 | 0.941 | 0.942 |
| scGPT32 线性 attention | 0.503 | 0.787 | 0.550 |
| scGPT32 CTA | 0.615 | 0.704 | 0.750 |
| scGPT32 HA | 0.554 | 0.627 | 0.508 |
| scGPT32 MixMIL | 0.887 | 0.911 | 0.692 |
| 原生 scGPT512 均值 | 0.923 | 0.941 | 0.808 |
| 原生 scGPT512 CTA | 0.605 | 0.822 | 0.450 |
| 原生 scGPT512 HA | 0.677 | 0.917 | 0.583 |

患者组织级 OOF 以病灶等权汇总。T/NK 为 16 位患者、28 个患者组织、12 个 PT–LM 对，Mono/Mac/DC 为 15、26、11，B/Plasma 为 15、23、8。每个方法的三个谱系合计 89 个袋，所有方法使用相同的样本、患者折和抽取细胞。其余分类指标与每个共同患者集的差值区间另附。[患者级完整指标](patient_metrics.csv) · [相对同表示均值头的配对区间](matched_mean_comparisons.csv) · [患者组织 logits](patient_site_logits.csv)

## 实现与验收

HA 和 CTA 以作者定义的 ReLU 投影、类型内 attention 或均值、类型间 attention 和线性分类头为基础。袋内缺失的 M2 子型不参与该袋的 softmax。单细胞的带符号贡献为类型内权重、类型间权重和该细胞经分类头的带符号投影之积。取相同权重及每袋均有两种亚型的算例，两个模型均与作者 PyG 实现的输出一致。32 维 HA/CTA 实际分别优化 1,155/1,122 个参数，原生 512 维为 16,515/16,482。[作者代码](https://github.com/minhchaudo/hier-mil/blob/2fb155840394961595755bede506c4b54c0942ae/model.py) · [方法论文](https://doi.org/10.1093/bioinformatics/btaf241)

MixMIL 使用作者发布的全协方差高斯变分后验及训练袋均值 GLMM 初始化，保留两组可学习先验尺度、截距和随机袋效应。其 32 维版本的 64 维联合高斯均值与对称协方差共有 2,144 个自由度；作者使用非三角 64×64 协方差因子，实际优化 4,160 个变分坐标，连同截距和两项尺度一共 4,163 个。训练阶段每次用 8 个后验样本估计 ELBO，并加入与普通头一致的患者配对项。OOF 采用后验均值参数代入的 logit，归一化尺度来自训练袋；其 sigmoid 对应代入参数的概率。固定参数、相同样本的袋随机效应与原作者张量函数最大绝对差为 2.38×10⁻⁷。[作者模型](https://github.com/AIH-SGML/mixmil/blob/16d92349ccf19b73c1dccd167934fcc3202b6be4/mixmil/model.py) · [MixMIL 论文](https://proceedings.mlr.press/v238/engelmann24a.html) · [数值对照](upstream_equivalence.json)

PCA32、scGPT32 和原生 scGPT512 分别实跑 45、75、30 个模型折，共 150 次。三个分支各用 CPU 2 线程，耗时分别为 31.92、33.26、18.51 秒，峰值 RSS 分别为 1.45、1.25、1.21 GiB。每个直接前作模型折抽取一个袋，核对单袋和批量推理一致。另对全部 148,880 行细胞贡献逐袋核验权重和为 1、带符号细胞贡献加偏置重构 logit。每个方法只有一个固定配置和一个种子。小开发集的折间模型各自拟合，bootstrap 区间条件于本次 OOF 预测。[设计与预算](design.md) · [PCA32 运行](PCA32/benchmark_runs.csv) · [scGPT32 运行](scGPT32/benchmark_runs.csv) · [原生 scGPT512 运行](scGPT512/benchmark_runs.csv)

## 重建

从项目根目录执行。依赖已通过 M4 的 scGPT 环境与 CPU 汇总环境安装，新增包版本另见 [extra_packages.txt](extra_packages.txt)。官方下载的软件保存在用户缓存，模型、投影和训练尺度位于 [models/m4_direct](../../../models/m4_direct/)。

```bash
cache="$HOME/.cache/crlm-m4-priors"
mkdir -p "$cache"
git clone https://github.com/minhchaudo/hier-mil "$cache/hier-mil"
git -C "$cache/hier-mil" checkout 2fb155840394961595755bede506c4b54c0942ae
git clone https://github.com/AIH-SGML/mixmil "$cache/mixmil"
git -C "$cache/mixmil" checkout 16d92349ccf19b73c1dccd167934fcc3202b6be4
uv pip install --python .venv-scgpt/bin/python --target "$cache/pylibs" --no-deps -r artifacts/m4/direct_prior/extra_packages.txt
export PYTHONPATH="$cache/pylibs" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
for rep in PCA32 scGPT32 scGPT512; do
    uv run --no-project --python .venv-scgpt/bin/python python -Bu scripts/benchmark_prior_mil.py --representation "$rep"
done
uv run --no-project --python .venv-scgpt/bin/python python -B artifacts/m4/direct_prior/check_upstream.py
PYTHONPATH=src uv run --no-project --python .venv-m4/bin/python python -Bm crlm.prior_mil_metrics
```

各分支目录保存全部 OOF 预测、细胞权重、带符号贡献、模型参数数及资源记录。源码、输入缓存、模型与输出 SHA 见 [清单](output_sha256.txt)。
