# 直接前作的固定比较

每个谱系沿用 M4 的患者五折、≥50 个细胞的袋和同一组每袋至多 256 个细胞。PCA 使用 M4 已保存的训练折基因、中心与 32 维分量。scGPT 使用同批冻结的 512 维细胞表示；这组直接比较在每个训练折内先投影到 32 维，保存中心与投影，HA、CTA、MixMIL 及相应线性 mean/attention 互用完全相同输入。这样 full-covariance MixMIL 的变分协方差因子规模固定在 64×64。主 M4 的原生 512 维线性头另列。折内表示逐维标准化，训练细胞提供所有变换参数。

HA、CTA 使用作者默认的 32 单元 ReLU 输入层、单一共享细胞层线性打分器、类型层打分器和线性输出头。HA 在同一亚型内学习细胞权重，CTA 对同一亚型等权。各自只对袋中出现的 M2 亚型归一化类型权重。使用作者 100 epoch、Adam 学习率 0.001、weight decay 0.001、dropout 0 的固定运行。实现用 PyTorch 的分组归约复写聚合，保留作者层级并对空亚型做显式 mask。[作者代码](https://github.com/minhchaudo/hier-mil/blob/2fb155840394961595755bede506c4b54c0942ae/model.py) · [论文](https://doi.org/10.1093/bioinformatics/btaf241)

MixMIL 保留作者全协方差高斯变分后验、两组各向同性高斯先验和可学习尺度、Bernoulli 袋似然及截距。使用官方变分后验源码；以训练袋均值的贝叶斯 GLMM 初始化，按作者代码的抽样归一化计算袋随机效应。预测使用变分均值代入，使用训练袋的均值/标准差固定预测尺度，并加入训练截距。每次 Adam 更新取 8 个后验样本，学习率 0.001，固定 500 步。作者软件可使用 2,000 步；本研究每个方法只试一个预先指定的配置，不查看外层标签调参。[官方后验](https://github.com/AIH-SGML/mixmil/blob/16d92349ccf19b73c1dccd167934fcc3202b6be4/mixmil/posterior.py) · [官方模型与归一化](https://github.com/AIH-SGML/mixmil/blob/16d92349ccf19b73c1dccd167934fcc3202b6be4/mixmil/model.py) · [论文](https://proceedings.mlr.press/v238/engelmann24a.html)

三种方法的 BCE 依患者、现有组织、病灶、技术样本层层等权，另加 0.1 倍患者内配对 softplus。MixMIL 的后验 KL 除以独立训练患者数。BCE 和配对损失均在同一组后验抽样上求期望。全训练袋批量、固定种子 20261121，固定 epoch 末读取；训练标签只进入训练过程。所有逐细胞预测贡献对应模型 logit 并另列截距。

另在完整 512 维 scGPT 缓存上运行 HA 与 CTA 并对照主 M4 的原生 scGPT mean。层级头的 32 单元隐藏层和优化设置相同。这一补充用于检查作者常用输入维度对固定模型的影响。
