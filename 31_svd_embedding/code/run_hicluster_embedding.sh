#!/bin/bash

# 运行 hicluster embedding 命令
# 使用生成的 cell_table.tsv 进行 SVD 嵌入分析
RESOLUTION=1000000
SCALE_FACTOR=100000
DIST=1000000
DIM=50
echo "开始运行 hicluster embedding..."

# 切换到 schicluster 环境并运行命令
/opt/homebrew/bin/micromamba run -n schicluster hicluster embedding \
    --cell_table_path "/Volumes/SumSung500/CSU/0_HiRES/31_svd_embedding/code/cell_table.tsv" \
    --output_dir "/Volumes/SumSung500/CSU/0_HiRES/31_svd_embedding/output/${DIM}d_${RESOLUTION}" \
    --dim ${DIM} \
    --dist ${DIST} \
    --resolution ${RESOLUTION} \
    --scale_factor ${SCALE_FACTOR} \
    --norm_sig \
    --save_raw \
    --cpu 20

echo "hicluster embedding 完成！"
echo "输出文件保存在: /Volumes/SumSung500/CSU/0_HiRES/31_svd_embedding/output"~