# 串行运行（低内存占用，但速度慢）
# python pair2h5ad_optimized.py \
#     --cell_table /path/to/your/cell_table.tsv \
#     --output /path/to/output.h5ad \
#     --chr_sizes /path/to/your/hg19.chrom.sizes \
#     --resolution 100000 \
#     --workers 1


# 并行运行（高内存占用，但速度快）
# 将 --workers 设置为一个大于1的数字。设置为 0 将自动使用你机器上所有的物理CPU核心，通常是最佳选择
# 使用所有可用的物理核心
python pair2h5ad_optimized.py \
    --cell_table /Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/20length_cell_table.tsv \
    --output /Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/output \
    --chr_sizes /Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/mm10_chrom_sizes.txt \
    --resolution 100000 \
    --workers 0