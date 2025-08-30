#!/bin/bash
# 1_run_add_metadata.sh - 运行添加元数据脚本
# 使用micromamba的3_schicluster_python38环境

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 打印带颜色的消息
print_message() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

print_message $BLUE "=================================================="
print_message $BLUE "为h5ad文件添加元数据"
print_message $BLUE "=================================================="

# 检查micromamba环境
print_message $BLUE "\n1. 检查micromamba环境..."
if ! command -v micromamba &> /dev/null; then
    print_message $RED "错误：未找到micromamba命令"
    exit 1
fi

# 激活schicluster环境
print_message $YELLOW "激活3_schicluster_python38环境..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate 3_schicluster_python38

# 检查Python环境
print_message $BLUE "\n2. 检查Python环境..."
python --version
which python

# 检查必要的Python包
print_message $BLUE "\n3. 检查必要的Python包..."
required_packages=("scanpy" "anndata" "pandas" "openpyxl")

for package in "${required_packages[@]}"; do
    if python -c "import ${package}" 2>/dev/null; then
        print_message $GREEN "✓ ${package} 已安装"
    else
        print_message $RED "✗ ${package} 未安装"
        print_message $YELLOW "正在安装 ${package}..."
        pip install "$package"
    fi
done

# 设置脚本路径和参数
SCRIPT_DIR="/home/duxuyan/Projects/schicluster_mac/5_plot_umap"
SCRIPT_FILE="${SCRIPT_DIR}/add_metadata_to_h5ad.py"
INPUT_DIR="/home/duxuyan/Projects/schicluster_mac/4_contact_decay_profile_2_h5ad/outputs"
METADATA_FILE="${SCRIPT_DIR}/GSE223917_HiRES_emb_metadata.xlsx"
OUTPUT_DIR="/home/duxuyan/Projects/schicluster_mac/4_contact_decay_profile_2_h5ad/outputs_with_metadata"

# 检查脚本文件
if [ ! -f "${SCRIPT_FILE}" ]; then
    print_message $RED "错误：找不到脚本文件: ${SCRIPT_FILE}"
    exit 1
fi

# 检查输入目录
if [ ! -d "${INPUT_DIR}" ]; then
    print_message $RED "错误：输入目录不存在: ${INPUT_DIR}"
    exit 1
fi

# 检查元数据文件
if [ ! -f "${METADATA_FILE}" ]; then
    print_message $RED "错误：元数据文件不存在: ${METADATA_FILE}"
    exit 1
fi

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 检查是否提供了stage参数
if [ ! -z "$1" ]; then
    STAGE_PARAM="--stage $1"
    CMD="python ${SCRIPT_FILE} --input_dir ${INPUT_DIR} --metadata_file ${METADATA_FILE} --output_dir ${OUTPUT_DIR} ${STAGE_PARAM}"
else
    CMD="python ${SCRIPT_FILE} --input_dir ${INPUT_DIR} --metadata_file ${METADATA_FILE} --output_dir ${OUTPUT_DIR}"
fi

# 显示将要执行的命令
print_message $BLUE "\n4. 执行命令:"
print_message $YELLOW "${CMD}"

# 记录开始时间
START_TIME=$(date +%s)
print_message $GREEN "\n开始时间: $(date)"

# 执行命令
eval $CMD

# 记录结束时间
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

print_message $GREEN "\n✅ 元数据添加完成！"
print_message $GREEN "结束时间: $(date)"
print_message $GREEN "总耗时: ${DURATION} 秒"

# 显示输出文件
print_message $BLUE "\n📁 输出文件:"
if [ -d "${OUTPUT_DIR}" ]; then
    ls -lh "${OUTPUT_DIR}"
else
    print_message $YELLOW "输出目录不存在: ${OUTPUT_DIR}"
fi

print_message $GREEN "\n🎉 元数据添加任务完成！"