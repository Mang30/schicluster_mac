#!/bin/bash
# Hi-C接触衰减特征UMAP分析运行脚本
# 作者：Claude Code Assistant  
# 日期：2025-08-29

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_message() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# 脚本目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"

print_message $BLUE "=================================================="
print_message $BLUE "Hi-C接触衰减特征UMAP分析工具"
print_message $BLUE "=================================================="

# 检查输入数据
DATA_FILE="${PROJECT_ROOT}/outputs/all_stages_decay_profiles.csv"
OUTPUT_DIR="${PROJECT_ROOT}/outputs/umap_analysis"

if [ ! -f "$DATA_FILE" ]; then
    print_message $RED "错误：未找到衰减曲线数据文件"
    print_message $YELLOW "请先运行接触衰减分析: ./run_contact_decay_analysis.sh"
    print_message $YELLOW "预期文件位置: $DATA_FILE"
    exit 1
fi

print_message $GREEN "✓ 找到输入数据: $DATA_FILE"

# 激活环境
print_message $BLUE "\n1. 激活schicluster环境..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate schicluster

# 检查和安装UMAP相关包
print_message $BLUE "\n2. 检查UMAP分析依赖..."
required_packages=("umap-learn" "scikit-learn")

for package in "${required_packages[@]}"; do
    if python -c "import ${package//-/_}" 2>/dev/null; then
        print_message $GREEN "✓ ${package} 已安装"
    else
        print_message $YELLOW "安装 ${package}..."
        pip install "$package"
    fi
done

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

print_message $BLUE "\n3. 开始UMAP分析..."
print_message $YELLOW "输入文件: $DATA_FILE"
print_message $YELLOW "输出目录: $OUTPUT_DIR"

# 设置参数
MIN_DISTANCE=100    # 最小距离 (kb)
MAX_DISTANCE=5000   # 最大距离 (kb) 
FEATURE_SELECTION="log_uniform"

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --min-distance)
            MIN_DISTANCE="$2"
            shift 2
            ;;
        --max-distance)
            MAX_DISTANCE="$2" 
            shift 2
            ;;
        --feature-selection)
            FEATURE_SELECTION="$2"
            shift 2
            ;;
        -h|--help)
            print_message $GREEN "用法: $0 [选项]"
            print_message $YELLOW "选项:"
            print_message $YELLOW "  --min-distance N        最小分析距离 (kb, 默认: 100)"
            print_message $YELLOW "  --max-distance N        最大分析距离 (kb, 默认: 5000)"
            print_message $YELLOW "  --feature-selection S   特征选择策略 (all/uniform/log_uniform, 默认: log_uniform)"
            print_message $YELLOW "  -h, --help              显示帮助"
            exit 0
            ;;
        *)
            print_message $RED "未知选项: $1"
            exit 1
            ;;
    esac
done

print_message $YELLOW "分析参数:"
print_message $YELLOW "  距离范围: ${MIN_DISTANCE}-${MAX_DISTANCE} kb"
print_message $YELLOW "  特征选择: ${FEATURE_SELECTION}"

# 运行UMAP分析
START_TIME=$(date +%s)
print_message $GREEN "\n开始时间: $(date)"

python "${SCRIPT_DIR}/scripts/create_umap_analysis.py" \
    --input "$DATA_FILE" \
    --output "$OUTPUT_DIR" \
    --min-distance "$MIN_DISTANCE" \
    --max-distance "$MAX_DISTANCE" \
    --feature-selection "$FEATURE_SELECTION"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

print_message $GREEN "\n✅ UMAP分析完成！"
print_message $GREEN "结束时间: $(date)"
print_message $GREEN "总耗时: ${DURATION} 秒"

# 显示结果
print_message $BLUE "\n📊 分析结果："
if [ -f "${OUTPUT_DIR}/cell_feature_matrix.csv" ]; then
    CELL_COUNT=$(tail -n +2 "${OUTPUT_DIR}/cell_feature_matrix.csv" | wc -l)
    FEATURE_COUNT=$(head -n 1 "${OUTPUT_DIR}/cell_feature_matrix.csv" | tr ',' '\n' | wc -l)
    print_message $YELLOW "细胞数量: $CELL_COUNT"
    print_message $YELLOW "特征维度: $((FEATURE_COUNT - 3))"  # 减去元数据列
fi

print_message $BLUE "\n📁 输出文件:"
if [ -d "$OUTPUT_DIR" ]; then
    find "$OUTPUT_DIR" -name "*.csv" -type f | while read -r file; do
        print_message $YELLOW "  数据: $(basename "$file")"
    done
    
    find "$OUTPUT_DIR" -name "*.png" -type f | while read -r file; do
        print_message $YELLOW "  图表: $(basename "$file")"
    done
fi

print_message $GREEN "\n🎉 UMAP分析完成！"
print_message $BLUE "\n💡 后续分析建议:"
print_message $YELLOW "1. 查看细胞特征矩阵: head ${OUTPUT_DIR}/cell_feature_matrix.csv"
print_message $YELLOW "2. 查看UMAP坐标: head ${OUTPUT_DIR}/umap_embedding.csv"  
print_message $YELLOW "3. 查看图表: open ${OUTPUT_DIR}/*.png"
print_message $YELLOW "4. 进一步聚类分析: 可使用特征矩阵进行聚类"