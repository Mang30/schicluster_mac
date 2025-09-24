#!/bin/bash

#SBATCH -o pair2h5ad.%j.out       ## Job's standard output and error file
#SBATCH -J pair2h5ad              ## Job name
#SBATCH -A pi_test                ## Account to charge
#SBATCH -p cpuQ                   ## Partition (queue) to use
#SBATCH -q cpuq                   ## QOS to use

# --- 关键资源请求 ---
# 您的 Python 脚本使用 multiprocessing，它在单个节点内利用多核心。
# 因此，我们只需要申请1个节点，但需要该节点上的大量CPU核心。
#SBATCH --nodes=1                 ## 申请1个计算节点
#SBATCH --ntasks-per-node=1         ## 在这个节点上只运行1个主任务 (即您的python脚本)
#SBATCH --cpus-per-task=48        ## 为这个任务申请48个CPU核心

# --- 打印作业信息 ---
echo "Job started on: $(date)"
echo "Running on node: $(hostname)"
echo "Allocated CPUs: $SLURM_CPUS_PER_TASK"
echo "Job ID: $SLURM_JOB_ID"
echo "--------------------"

# --- 准备运行环境 ---
# !! 重要 !!: 在Slurm脚本中激活conda/mamba环境
# 必须先初始化shell，然后再激活环境
source /path/to/your/micromamba/etc/profile.d/micromamba.sh # 替换为您的 micromamba 初始化路径
micromamba activate schicluster

# 验证Python环境是否正确
echo "Using python executable at: $(which python)"
echo "--------------------"

# --- 定义文件路径 (推荐使用变量，方便修改) ---
CELL_TABLE="/Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/20length_cell_table.tsv"
OUTPUT_H5AD="/Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/output/my_output.h5ad" # 确保输出是文件而不是目录
CHR_SIZES="/Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/mm10_chrom_sizes.txt"

# --- 运行您的Python脚本 ---
# 我们将Slurm分配的核心数 ($SLURM_CPUS_PER_TASK) 直接传递给脚本的 --workers 参数
# 这样脚本就能完美利用所有申请到的资源
python pair2h5ad_optimized.py \
    --cell_table ${CELL_TABLE} \
    --output ${OUTPUT_H5AD} \
    --chr_sizes ${CHR_SIZES} \
    --resolution 100000 \
    --workers $SLURM_CPUS_PER_TASK

echo "--------------------"
echo "Job finished on: $(date)"



#
# sbatch -A pi_test -p cpuQ -q cpuq run_pair2h5ad.sh