#!/bin/bash

# run_stage_E75_umap.sh - 运行Stage E75 UMAP绘图脚本

echo "=== Stage E75 UMAP 绘图脚本 ==="
echo "开始时间: $(date)"

# 切换到脚本目录
cd "$(dirname "$0")"

# 激活conda环境
echo "激活schicluster环境..."
source /Users/wuhaoliu/mamba/etc/profile.d/mamba.sh
micromamba activate schicluster

# 检查Python环境
echo "Python路径: $(which python)"
echo "Python版本: $(python --version)"

# 检查必要的库
echo "检查必要的库..."
python -c "
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
print('所有必要的库都已安装')
print(f'scanpy版本: {sc.__version__}')
print(f'pandas版本: {pd.__version__}')
print(f'numpy版本: {np.__version__}')
"

if [ $? -ne 0 ]; then
    echo "错误: 缺少必要的Python库"
    exit 1
fi

# 检查输入文件
input_file="/Volumes/SumSung500/CSU/0_HiRES/4_contact_decay_profile_2_h5ad/json2h5ad/output/stage_E75.h5ad"
if [ ! -f "$input_file" ]; then
    echo "错误: 输入文件不存在: $input_file"
    echo "请确保已经运行了JSON到H5AD的转换脚本"
    exit 1
fi

echo "输入文件检查通过: $input_file"

# 运行Python脚本
echo "开始运行UMAP绘图脚本..."
python plot_stage_E75_umap.py

# 检查执行结果
if [ $? -eq 0 ]; then
    echo "=== UMAP绘图完成 ==="
    echo "结束时间: $(date)"
    
    # 显示输出文件
    output_dir="./output"
    if [ -d "$output_dir" ]; then
        echo "输出文件列表:"
        ls -la "$output_dir"/*.png 2>/dev/null || echo "没有找到PNG文件"
        ls -la "$output_dir"/*.txt 2>/dev/null || echo "没有找到文本文件"
    fi
else
    echo "错误: UMAP绘图脚本执行失败"
    exit 1
fi
