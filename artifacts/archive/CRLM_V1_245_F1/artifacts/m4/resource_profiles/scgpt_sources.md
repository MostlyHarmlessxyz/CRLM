# scGPT whole-human 编码环境

使用作者发布的 whole-human 编码器架构和冻结参数。预训练说明为约 3,300 万正常人类细胞；逐细胞训练清单未取得，预训练与研究样本的重叠程度未知。源码固定为 `cebd6fae655b9c585a4807daa3ac31bb764f06b4`。[作者源码及模型目录](https://github.com/bowang-lab/scGPT/tree/cebd6fae655b9c585a4807daa3ac31bb764f06b4)

本机从 wanglab/scGPT-human 修订 `a24c237737a40f3720f75abb555489e9fe753be6` 取得 args.json、vocab.json、best_model.pt。三份文件的实际 SHA-256 见 [profile](scgpt_profile.json)，与独立研究给出的 whole-human 发布文件指纹逐件相符。[固定权重版本](https://huggingface.co/wanglab/scGPT-human/tree/a24c237737a40f3720f75abb555489e9fe753be6) · [核对来源](https://doi.org/10.1038/s41598-026-77853-3)

项目入口加载固定源码中的 TransformerModel，并按作者非 flash 后端的参数名转换核对每个实际参与 CLS 编码的张量。输入按库读取作者交付的稀疏整数，先匹配词表，再对每个细胞做 51 档分箱和固定随机流的 token 抽样，最多 1,199 个表达 token 加一个 CLS。输出为 512 维 L2 归一化 CLS 向量。每个缓存保存细胞 ID、测到的非零基因数、词表覆盖数、实际 token 数、输入与实现版本、编码耗时。

在项目根目录可按下面的命令准备与实测相同的输入。此处使用独立的 PyTorch XPU 环境，CPU 参照也由同一份 torch 运行。

```bash
uv venv --python 3.12.13 .venv-scgpt
uv pip install --python .venv-scgpt/bin/python --index-url https://download.pytorch.org/whl/xpu 'torch==2.15.0+xpu'
uv pip install --python .venv-scgpt/bin/python -r artifacts/m1/python-packages.txt 'tqdm==4.70.1'
model_dir=models/scGPT-whole-human
mkdir -p "$model_dir"
revision=a24c237737a40f3720f75abb555489e9fe753be6
for file in args.json vocab.json best_model.pt; do
  curl -fLsS "https://huggingface.co/wanglab/scGPT-human/resolve/$revision/$file" -o "$model_dir/$file" || exit 1
done
commit=cebd6fae655b9c585a4807daa3ac31bb764f06b4
curl -fLsS "https://codeload.github.com/bowang-lab/scGPT/tar.gz/$commit" -o "$model_dir/scGPT_source.tar.gz"
tar -xzf "$model_dir/scGPT_source.tar.gz" -C "$model_dir"
uv run --no-project --python .venv-scgpt/bin/python python -Bu scripts/encode_scgpt.py --profile --batch 4
```
