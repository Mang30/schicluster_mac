#!/bin/bash
# Hi-C接触衰减曲线分析运行脚本 - stage_organized版本
# 处理hires_data_processing/outputs/stage_organized中的数据
# 使用micromamba的schicluster环境

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

STAGE_ORGANIZED_DIR="$(cd "${PROJECT_ROOT}/../hires_data_processing/outputs/stage_organized" 2>/dev/null && pwd || echo "${PROJECT_ROOT}/../hires_data_processing/outputs/stage_organized")"

BATCH_PROCESS_SCRIPT="${SCRIPT_DIR}/batch_process_stage_organized.py"

OUTPUT_DIR="${PROJECT_ROOT}/outputs"

# --- 默认参数 ---
STAGES="" # 空表示处理所有stage
MAX_WORKERS=6 # 默认使用所有可用的CPU核心
MAX_FILES="" # 默认不限制文件数

# --- 解析命令行参数 ---
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--stages)
            STAGES="$2"; shift 2 ;;
        -w|--max-workers)
            MAX_WORKERS="$2"; shift 2 ;;
        -n|--max-files)
            MAX_FILES="$2"; shift 2 ;;
        -h|--help)
            print_message $GREEN "用法: $0 [选项]"
            print_message $YELLOW "选项:"
            print_message $YELLOW "  -s, --stages STAGES      指定处理的阶段，用逗号分隔 (默认: 全部)"
            print_message $YELLOW "                           示例: E75,E85,EX05 或 all"
            print_message $YELLOW "  -w, --max-workers N      并行worker数量 (默认: all available cores)"
            print_message $YELLOW "  -n, --max-files N        每个阶段最大处理文件数 (默认: 不限制)"
            print_message $YELLOW "  -h, --help               显示帮助信息"
            print_message $GREEN "\n示例:"
            print_message $YELLOW "  # 处理所有stage的数据"
            print_message $YELLOW "  $0"
            print_message $YELLOW "  # 处理指定的stage"
            print_message $YELLOW "  $0 -s E75,E85,EX05"
            print_message $YELLOW "  # 处理单个stage，使用16个worker"
            print_message $YELLOW "  $0 -s E75 -w 16"
            print_message $YELLOW "  # 每个stage限制处理100个文件"
            print_message $YELLOW "  $0 -s E75,E85 -n 100"
            exit 0 ;;
        *)
            print_message $RED "未知选项: $1"; exit 1 ;;
    esac
done

# --- 脚本主体 ---
print_message $BLUE "=================================================="
print_message $BLUE "Hi-C接触衰减曲线分析工具 - Stage Organized版本"
print_message $BLUE "=================================================="

# 激活环境并检查依赖
print_message $BLUE "\n1. 准备环境..."
eval "$(micromamba shell hook --shell bash)"
#micromamba activate 3_schicluster_python38
micromamba activate schicluster
#print_message $GREEN "✓ 环境 '3_schicluster_python38' 已激活"
print_message $GREEN "✓ 环境 'schicluster' 已激活"

# 显示处理参数
print_message $BLUE "\n2. 分析参数配置:"
print_message $YELLOW "  Stage组织数据目录: ${STAGE_ORGANIZED_DIR}"
if [ -n "$STAGES" ] && [ "$STAGES" != "all" ]; then
    print_message $YELLOW "  处理阶段: ${STAGES}"
else
    print_message $YELLOW "  处理阶段: 全部可用阶段"
fi
print_message $YELLOW "  并行 Worker 数: ${MAX_WORKERS}"
if [ -n "$MAX_FILES" ]; then
    print_message $YELLOW "  每个阶段最大文件数: ${MAX_FILES}"
else
    print_message $YELLOW "  文件数量: 不限制（全量处理）"
fi

# 检查输入目录
print_message $BLUE "\n3. 检查输入数据路径..."
print_message $YELLOW "  预期路径: ${STAGE_ORGANIZED_DIR}"
if [ ! -d "${STAGE_ORGANIZED_DIR}" ]; then
    print_message $RED "错误：Stage组织数据目录不存在: ${STAGE_ORGANIZED_DIR}"
    print_message $YELLOW "提示：请先运行 organize_cool_files_by_stage.py 脚本组织数据"
    exit 1
fi
print_message $GREEN "✓ Stage组织数据目录存在"

# 检查有多少个stage目录
AVAILABLE_STAGES=($(ls -d ${STAGE_ORGANIZED_DIR}/*/  2>/dev/null | xargs -n 1 basename || echo ""))
if [ ${#AVAILABLE_STAGES[@]} -eq 0 ]; then
    print_message $RED "错误：在 ${STAGE_ORGANIZED_DIR} 中未找到任何stage目录"
    exit 1
fi
print_message $GREEN "✓ 发现 ${#AVAILABLE_STAGES[@]} 个stage目录: ${AVAILABLE_STAGES[*]}"

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"

# 构建命令
CMD="python '${BATCH_PROCESS_SCRIPT}' \
    --input '${STAGE_ORGANIZED_DIR}' \
    --output '${OUTPUT_DIR}' \
    --max-workers ${MAX_WORKERS}"

if [ -n "$STAGES" ] && [ "$STAGES" != "all" ]; then
    CMD="${CMD} --stages '${STAGES}'"
fi

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
print_message $BLUE "\n📁 请查看输出目录: ${OUTPUT_DIR}"
print_message $BLUE "   每个stage的结果保存在对应的子目录中"