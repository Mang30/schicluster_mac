#!/bin/bash
# Hi-C接触衰减曲线分析运行脚本 - 优化版
# 使用micromamba的schicluster环境
#
# 这个脚本现在更加灵活，允许通过命令行参数指定
# 要处理的阶段(stage)和分辨率(resolution)。

set -e # 遇到错误立即退出

# --- 颜色定义 ---
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

# --- 路径定义 ---
# 获取脚本的绝对路径，确保无论从哪个目录执行都能正确找到文件
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_ROOT="$(cd "${PROJECT_ROOT}/../hires_data_processing/outputs" 2>/dev/null && pwd || echo "${PROJECT_ROOT}/../hires_data_processing/outputs")"
BATCH_PROCESS_SCRIPT="${SCRIPT_DIR}/batch_process_stages.py"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

# --- 默认参数 ---
STAGE="EX05"
RESOLUTION="100K"
MAX_WORKERS=$(nproc --all) # 默认使用所有可用的CPU核心
MAX_FILES="" # 默认不限制文件数

# --- 解析命令行参数 ---
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--stage)
            STAGE="$2"; shift 2 ;;
        -r|--resolution)
            RESOLUTION="$2"; shift 2 ;;
        -w|--max-workers)
            MAX_WORKERS="$2"; shift 2 ;;
        -n|--max-files)
            MAX_FILES="$2"; shift 2 ;;
        -h|--help)
            print_message $GREEN "用法: $0 [选项]"
            print_message $YELLOW "选项:"
            print_message $YELLOW "  -s, --stage STAGE        指定处理的阶段 (默认: ${STAGE})"
            print_message $YELLOW "  -r, --resolution RES     指定分辨率目录 (默认: ${RESOLUTION})"
            print_message $YELLOW "  -w, --max-workers N      并行worker数量 (默认: all available cores)"
            print_message $YELLOW "  -n, --max-files N        每个阶段最大处理文件数 (默认: 不限制)"
            print_message $YELLOW "  -h, --help               显示帮助信息"
            print_message $GREEN "\n示例:"
            print_message $YELLOW "  # 处理 E75 阶段的 200K 分辨率数据，使用 16 个 worker"
            print_message $YELLOW "  $0 -s E75 -r 200K -w 16"
            exit 0 ;;
        *)
            print_message $RED "未知选项: $1"; exit 1 ;;
    esac
done

# --- 脚本主体 ---
print_message $BLUE "=================================================="
print_message $BLUE "Hi-C接触衰减曲线分析工具"
print_message $BLUE "=================================================="

# 激活环境并检查依赖
print_message $BLUE "\n1. 准备环境..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate 3_schicluster_python38
print_message $GREEN "✓ 环境 '3_schicluster_python38' 已激活"

# 显示处理参数
print_message $BLUE "\n2. 分析参数配置:"
print_message $YELLOW "  数据根目录: ${DATA_ROOT}"
print_message $YELLOW "  处理阶段: ${STAGE}"
print_message $YELLOW "  处理分辨率: ${RESOLUTION}"
print_message $YELLOW "  并行 Worker 数: ${MAX_WORKERS}"
if [ -n "$MAX_FILES" ]; then
    print_message $YELLOW "  最大文件数: ${MAX_FILES}"
else
    print_message $YELLOW "  文件数量: 不限制（全量处理）"
fi

# 检查输入目录
STAGE_DATA_PATH="${DATA_ROOT}/${STAGE}/impute/${RESOLUTION}"
print_message $BLUE "\n3. 检查输入数据路径..."
print_message $YELLOW "  预期路径: ${STAGE_DATA_PATH}"
if [ ! -d "${STAGE_DATA_PATH}" ]; then
    print_message $RED "错误：输入数据目录不存在: ${STAGE_DATA_PATH}"
    exit 1
fi
print_message $GREEN "✓ 输入数据目录存在"

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 构建命令
CMD="python '${BATCH_PROCESS_SCRIPT}' \
    --input '${DATA_ROOT}' \
    --output '${OUTPUT_DIR}' \
    --stage '${STAGE}' \
    --resolution '${RESOLUTION}' \
    --max-workers ${MAX_WORKERS}"

if [ -n "$MAX_FILES" ]; then
    CMD="${CMD} --max-files ${MAX_FILES}"
fi

# 执行命令
print_message $BLUE "\n4. 开始执行分析..."
print_message $YELLOW "执行命令: ${CMD}"

START_TIME=$(date +%s)
eval $CMD
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# 结果展示
print_message $GREEN "\n✅ 分析完成！"
print_message $GREEN "总耗时: ${DURATION} 秒"
print_message $BLUE "\n📁 请查看输出目录: ${OUTPUT_DIR}/stage_${STAGE}_${RESOLUTION}"
