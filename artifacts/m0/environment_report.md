# 本机 CPU 环境

2026-11-21 11:43 UTC 的 [系统探测](environment_probe.json) 记录了 Arch Linux、内核 7.3.5-arch1-1、i5-13400 和 16 个逻辑 CPU。系统可见内存 62.56 GiB，当时可用 54.27 GiB。项目位于本地 NVMe 的 Btrfs，探测时可用空间 480.22 GiB。M1 按库串行处理，沿用配置中的 80 GB 空闲空间规划值。

系统 Python 为 3.14.10，NumPy 2.5.5、pandas 2.3.3，未安装 SciPy、AnnData、Scanpy、torch 或 scGPT，PATH 中无 R。uv 0.13.9 在用户目录建立独立 starter 环境，解释器 3.12.13，历史自检使用 3.12.14。当前依赖精确版本见 [环境清单](baseline_checks/python-packages.txt)。

`lspci -nnk` 显示 A750 和 UHD 730 均由 i915 驱动。加载已安装 oneAPI 的 compiler、umf、tbb 运行库后，`sycl-ls` 在 Level Zero 和 OpenCL 均枚举到独立 A750。torch XPU、scGPT 数值与性能待 M4 小批实测。M0 使用本地 CPU。
