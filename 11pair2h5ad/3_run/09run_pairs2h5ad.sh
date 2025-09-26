#!/bin/bash

#SBATCH -o ../logs/pairs2h5ad/9chr_pairs2h5ad.%j.out       ## Job's standard output and error file
#SBATCH -J 9chr_pairs2h5ad              ## Job name
#SBATCH -A free                   ## Account to charge
#SBATCH -p freecpuQ               ## Partition (queue) to use
#SBATCH -q normal                 ## QOS to use

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
# 直接定义Python解释器的绝对路径
PYTHON_EXEC="/public/home/hpc254701055/micromamba/envs/3_schicluster_python38/bin/python"

echo "Using python executable directly at: ${PYTHON_EXEC}"
echo "--------------------"

# --- 定义文件路径 (推荐使用变量，方便修改) ---
INPUT_DIR="/public/home/hpc254701055/Downloads/11pair2h5ad/output/split_pairs"
OUTPUT_DIR="/public/home/hpc254701055/Downloads/11pair2h5ad/output/chr_h5ad_${RESOLUTION}/" # 确保输出是文件而不是目录
CHR_SIZES="/public/home/hpc254701055/Downloads/11pair2h5ad/mm10_chrom_sizes_with_chrY.txt"
WORKERS=48
RESOLUTION=1000000

# --- 运行您的Python脚本 ---
# 我们将Slurm分配的核心数 ($SLURM_CPUS_PER_TASK) 直接传递给脚本的 --workers 参数
# 这样脚本就能完美利用所有申请到的资源
${PYTHON_EXEC}  ../chrData_to_h5ad_code/oneByone_convert_intra_pairs_to_h5ad.py\
    -i ${INPUT_DIR} \
    -s ${CHR_SIZES} \
    -c chr9 \
    -o ${OUTPUT_DIR} \
    -w ${WORKERS}\
    -r ${RESOLUTION}

echo "--------------------"
echo "Job finished on: $(date)"


