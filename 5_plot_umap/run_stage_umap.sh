#!/bin/bash

# run_stage_umap.sh - 通用的Stage UMAP绘图脚本

# 使用方法:
# ./run_stage_umap.sh stage_E75
# ./run_stage_umap.sh stage_E70
# 等等...

if [ $# -eq 0 ]; then
    echo "使用方法: $0 <stage_name> [h5ad_file] [output_dir]"
    echo "例如: $0 stage_E75"
    echo "     $0 stage_E70 /path/to/stage_E70.h5ad /path/to/output"
    exit 1
fi

STAGE_NAME=$1
H5AD_FILE=$2
OUTPUT_DIR=$3

echo "=== $STAGE_NAME UMAP 绘图脚本 ==="
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
" > /dev/null 2>&1

if [ $? -ne 0 ]; then
    echo "错误: 缺少必要的Python库"
    exit 1
fi

# 检查输入文件（如果没有指定，使用默认路径）
if [ -z "$H5AD_FILE" ]; then
    default_input_file="/Volumes/SumSung500/CSU/0_HiRES/4_contact_decay_profile_2_h5ad/json2h5ad/output/${STAGE_NAME}.h5ad"
    if [ ! -f "$default_input_file" ]; then
        echo "错误: 默认输入文件不存在: $default_input_file"
        echo "请确保已经运行了JSON到H5AD的转换脚本，或者指定正确的H5AD文件路径"
        exit 1
    fi
    echo "使用默认输入文件: $default_input_file"
else
    if [ ! -f "$H5AD_FILE" ]; then
        echo "错误: 指定的输入文件不存在: $H5AD_FILE"
        exit 1
    fi
    echo "使用指定输入文件: $H5AD_FILE"
fi

# 构建命令
CMD="python plot_stage_umap.py --stage $STAGE_NAME"

if [ ! -z "$H5AD_FILE" ]; then
    CMD="$CMD --h5ad_file $H5AD_FILE"
fi

if [ ! -z "$OUTPUT_DIR" ]; then
    CMD="$CMD --output_dir $OUTPUT_DIR"
fi

# 运行Python脚本
echo "开始运行UMAP绘图脚本..."
echo "执行命令: $CMD"
eval $CMD

# 检查执行结果
if [ $? -eq 0 ]; then
    echo "=== $STAGE_NAME UMAP绘图完成 ==="
    echo "结束时间: $(date)"
    
    # 显示输出文件
    if [ -z "$OUTPUT_DIR" ]; then
        output_path="./output/$STAGE_NAME"
    else
        output_path="$OUTPUT_DIR/$STAGE_NAME"
    fi
    
    if [ -d "$output_path" ]; then
        echo "输出文件列表 ($output_path):"
        ls -la "$output_path"/*.png 2>/dev/null || echo "没有找到PNG文件"
        ls -la "$output_path"/*.txt 2>/dev/null || echo "没有找到文本文件"
    fi
else
    echo "错误: $STAGE_NAME UMAP绘图脚本执行失败"
    exit 1
fi
