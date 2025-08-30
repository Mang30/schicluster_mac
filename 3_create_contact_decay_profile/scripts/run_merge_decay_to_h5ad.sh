#!/bin/bash
# 将E75阶段单细胞衰减曲线数据合并到h5ad文件并进行UMAP可视化
# 作者：Claude Code Assistant
# 日期：2025-08-29

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

# 脚本目录和项目根目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

print_message $BLUE "=================================================="
print_message $BLUE "E75阶段衰减曲线数据合并到h5ad文件"
print_message $BLUE "=================================================="
print_message $YELLOW "脚本目录: ${SCRIPT_DIR}"
print_message $YELLOW "项目根目录: ${PROJECT_ROOT}"

# 检查micromamba环境
print_message $BLUE "\n1. 检查micromamba环境..."
if ! command -v micromamba &> /dev/null; then
    print_message $RED "错误：未找到micromamba命令"
    exit 1
fi

# 激活schicluster环境
print_message $YELLOW "激活schicluster环境..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate 3_schicluster_python38

# 检查Python环境
print_message $BLUE "\n2. 检查Python环境..."
python --version
which python

# 检查必要的Python包
print_message $BLUE "\n3. 检查必要的Python包..."
required_packages=("scanpy" "anndata" "pandas" "numpy" "matplotlib")

for package in "${required_packages[@]}"; do
    if python -c "import ${package}" 2>/dev/null; then
        print_message $GREEN "✓ ${package} 已安装"
    else
        print_message $RED "✗ ${package} 未安装"
        print_message $YELLOW "正在安装 ${package}..."
        pip install "$package"
    fi
done

# 设置输入输出目录
INPUT_DIR="${PROJECT_ROOT}/outputs/stage_E75"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

print_message $BLUE "\n4. 处理参数..."
print_message $YELLOW "输入目录: ${INPUT_DIR}"
print_message $YELLOW "输出目录: ${OUTPUT_DIR}"

# 检查输入目录
if [ ! -d "${INPUT_DIR}" ]; then
    print_message $RED "错误：输入目录不存在: ${INPUT_DIR}"
    exit 1
fi

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 构建命令
CMD="python ${SCRIPT_DIR}/merge_decay_to_h5ad.py"
CMD="${CMD} --input_dir '${INPUT_DIR}'"
CMD="${CMD} --output_dir '${OUTPUT_DIR}'"

# 显示将要执行的命令
print_message $BLUE "\n5. 执行命令:"
print_message $YELLOW "${CMD}"

# 记录开始时间
START_TIME=$(date +%s)
print_message $GREEN "\n开始时间: $(date)"

# 执行命令
eval $CMD

# 记录结束时间
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

print_message $GREEN "\n✅ E75衰减曲线数据合并完成！"
print_message $GREEN "结束时间: $(date)"
print_message $GREEN "总耗时: ${DURATION} 秒"

# 显示结果文件
print_message $BLUE "\n📁 输出文件:"
if [ -f "${OUTPUT_DIR}/E75_decay_profiles.h5ad" ]; then
    print_message $GREEN "✓ h5ad文件: ${OUTPUT_DIR}/E75_decay_profiles.h5ad"
    ls -lh "${OUTPUT_DIR}/E75_decay_profiles.h5ad"
fi

if [ -f "${OUTPUT_DIR}/E75_decay_summary.txt" ]; then
    print_message $GREEN "✓ 摘要文件: ${OUTPUT_DIR}/E75_decay_summary.txt"
    cat "${OUTPUT_DIR}/E75_decay_summary.txt"
fi

print_message $GREEN "\n🎉 E75阶段衰减曲线数据合并完成！"
print_message $YELLOW "生成的h5ad文件可用于UMAP可视化和其他下游分析。"