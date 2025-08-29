#!/bin/bash
# Hi-C接触衰减曲线分析运行脚本
# 使用micromamba的schicluster环境
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
PROJECT_ROOT="${SCRIPT_DIR}"
# 数据根目录：从3_create_contact_decay_profile到hires_data_processing/outputs
DATA_ROOT="${SCRIPT_DIR}/../hires_data_processing/outputs"

print_message $BLUE "=================================================="
print_message $BLUE "Hi-C接触衰减曲线分析工具"
print_message $BLUE "=================================================="
print_message $YELLOW "项目根目录: ${PROJECT_ROOT}"
print_message $YELLOW "数据根目录: ${DATA_ROOT}"
print_message $YELLOW "脚本目录: ${SCRIPT_DIR}"

# 检查数据目录
if [ ! -d "${DATA_ROOT}" ]; then
    print_message $RED "错误：数据目录不存在: ${DATA_ROOT}"
    exit 1
fi

# 检查micromamba环境
print_message $BLUE "\n1. 检查micromamba环境..."
if ! command -v micromamba &> /dev/null; then
    print_message $RED "错误：未找到micromamba命令"
    exit 1
fi

# 激活schicluster环境
print_message $YELLOW "激活schicluster环境..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate schicluster

# 检查Python环境
print_message $BLUE "\n2. 检查Python环境..."
python --version
which python

# 检查必要的Python包
print_message $BLUE "\n3. 检查必要的Python包..."
required_packages=("numpy" "pandas" "matplotlib" "seaborn" "cooler" "cooltools")

for package in "${required_packages[@]}"; do
    if python -c "import ${package}" 2>/dev/null; then
        print_message $GREEN "✓ ${package} 已安装"
    else
        print_message $RED "✗ ${package} 未安装"
        print_message $YELLOW "正在安装 ${package}..."
        
        # 尝试使用不同的安装方法
        if [ "$package" == "cooler" ] || [ "$package" == "cooltools" ]; then
            micromamba install -c conda-forge -c bioconda "$package" -y
        else
            pip install "$package"
        fi
    fi
done

# 创建输出目录
OUTPUT_DIR="${PROJECT_ROOT}/outputs"
mkdir -p "${OUTPUT_DIR}"

print_message $BLUE "\n4. 开始批量处理..."
print_message $YELLOW "输出目录: ${OUTPUT_DIR}"

# 设置默认参数
MAX_FILES_PER_STAGE=5  # 每个阶段最多处理5个文件（测试用）
MAX_WORKERS=2          # 并行处理worker数量
MAX_DISTANCE=1000      # 最大分析距离（bin数）
SPECIFIC_STAGE=""      # 指定处理的阶段

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--stage)
            SPECIFIC_STAGE="$2"
            shift 2
            ;;
        -n|--max-files)
            MAX_FILES_PER_STAGE="$2"
            shift 2
            ;;
        -w|--max-workers)
            MAX_WORKERS="$2"
            shift 2
            ;;
        -d|--max-distance)
            MAX_DISTANCE="$2"
            shift 2
            ;;
        --full)
            MAX_FILES_PER_STAGE=""  # 处理所有文件
            shift
            ;;
        --test)
            MAX_FILES_PER_STAGE=2   # 测试模式，每个阶段只处理2个文件
            shift
            ;;
        --list-stages)
            print_message $BLUE "可用的发育阶段:"
            if [ -d "${DATA_ROOT}" ]; then
                find "${DATA_ROOT}" -type d -name "E*" -o -name "EX*" | sort | while read -r stage_dir; do
                    stage_name=$(basename "$stage_dir")
                    print_message $YELLOW "  ${stage_name}"
                done
            else
                print_message $RED "数据目录不存在: ${DATA_ROOT}"
            fi
            exit 0
            ;;
        -h|--help)
            print_message $GREEN "用法: $0 [选项]"
            print_message $YELLOW "选项:"
            print_message $YELLOW "  -s, --stage STAGE      指定处理的阶段 (如: E70, E75)"
            print_message $YELLOW "  -n, --max-files N      每个阶段最大处理文件数 (默认: 5)"
            print_message $YELLOW "  -w, --max-workers N    并行worker数量 (默认: 2)"
            print_message $YELLOW "  -d, --max-distance N   最大分析距离 (默认: 1000)"
            print_message $YELLOW "  --full                 处理所有文件"
            print_message $YELLOW "  --test                 测试模式 (每阶段2个文件)"
            print_message $YELLOW "  --list-stages          列出所有可用的发育阶段"
            print_message $YELLOW "  -h, --help             显示帮助信息"
            print_message $GREEN "\n示例:"
            print_message $YELLOW "  $0 --stage E70         # 只处理E70阶段"
            print_message $YELLOW "  $0 --stage E75 --test  # 测试模式处理E75阶段"
            print_message $YELLOW "  $0 --list-stages       # 列出所有可用阶段"
            exit 0
            ;;
        *)
            print_message $RED "未知选项: $1"
            exit 1
            ;;
    esac
done

# 显示处理参数
print_message $YELLOW "处理参数:"
if [ -n "$SPECIFIC_STAGE" ]; then
    print_message $YELLOW "  指定阶段: ${SPECIFIC_STAGE}"
else
    print_message $YELLOW "  处理所有阶段"
fi
if [ -n "$MAX_FILES_PER_STAGE" ]; then
    print_message $YELLOW "  每阶段最大文件数: ${MAX_FILES_PER_STAGE}"
else
    print_message $YELLOW "  每阶段最大文件数: 不限制"
fi
print_message $YELLOW "  并行worker数: ${MAX_WORKERS}"
print_message $YELLOW "  最大分析距离: ${MAX_DISTANCE} bins"

# 构建命令
CMD="python ${SCRIPT_DIR}/scripts/batch_process_stages.py"
CMD="${CMD} --input '${DATA_ROOT}'"
CMD="${CMD} --output '${OUTPUT_DIR}'"
CMD="${CMD} --max-workers ${MAX_WORKERS}"

if [ -n "$SPECIFIC_STAGE" ]; then
    CMD="${CMD} --stage ${SPECIFIC_STAGE}"
fi

if [ -n "$MAX_FILES_PER_STAGE" ]; then
    CMD="${CMD} --max-files ${MAX_FILES_PER_STAGE}"
fi

if [ -n "$MAX_DISTANCE" ]; then
    CMD="${CMD} --max-distance ${MAX_DISTANCE}"
fi

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

print_message $GREEN "\n✅ 批量处理完成！"
print_message $GREEN "结束时间: $(date)"
print_message $GREEN "总耗时: ${DURATION} 秒"

# 显示结果摘要
if [ -f "${OUTPUT_DIR}/batch_processing_report.json" ]; then
    print_message $BLUE "\n📊 处理结果摘要:"
    
    # 提取关键信息
    TOTAL_STAGES=$(python -c "
import json
try:
    with open('${OUTPUT_DIR}/batch_processing_report.json', 'r') as f:
        data = json.load(f)
    print(f\"发育阶段数: {data.get('stages_processed', 0)}\")
    print(f\"成功分析: {data.get('total_successful', 0)}/{data.get('total_files', 0)}\")
    print(f\"成功率: {data.get('overall_success_rate', 0)*100:.1f}%\")
    print(f\"数据点数: {data.get('combined_data_points', 0)}\")
    
    # 显示各阶段结果
    stage_results = data.get('stage_results', {})
    for stage, result in stage_results.items():
        success_rate = result.get('successful', 0) / max(result.get('total_files', 1), 1) * 100
        print(f\"  {stage}: {result.get('successful', 0)}/{result.get('total_files', 0)} ({success_rate:.1f}%)\")
        
except Exception as e:
    print(f'无法读取报告文件: {e}')
" 2>/dev/null)
    
    echo -e "$TOTAL_STAGES"
fi

# 显示输出文件
print_message $BLUE "\n📁 输出文件:"
print_message $YELLOW "主要输出目录: ${OUTPUT_DIR}"

if [ -d "${OUTPUT_DIR}" ]; then
    find "${OUTPUT_DIR}" -name "*.png" -type f | head -5 | while read -r file; do
        print_message $YELLOW "  图表: $(basename "$file")"
    done
    
    find "${OUTPUT_DIR}" -name "*_summary.json" -type f | head -3 | while read -r file; do
        print_message $YELLOW "  摘要: $(basename "$file")"
    done
fi

print_message $GREEN "\n🎉 分析完成！请查看输出目录中的结果文件。"

# 提供后续操作建议
print_message $BLUE "\n💡 后续操作建议:"
print_message $YELLOW "1. 查看批量处理报告: cat ${OUTPUT_DIR}/batch_processing_report.json"
print_message $YELLOW "2. 查看对比图表: ls ${OUTPUT_DIR}/*.png"
print_message $YELLOW "3. 查看合并数据: head ${OUTPUT_DIR}/all_stages_decay_profiles.csv"
print_message $YELLOW "4. 分析单个样本: python ${SCRIPT_DIR}/src/contact_decay_analyzer.py -i <cool_file> -o <output_dir>"