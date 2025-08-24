#!/bin/bash

# 🚀 HiRES Hi-C 数据处理一键启动脚本
# 作者: Claude AI Assistant
# 版本: v1.0
# 日期: 2025-08-24

set -e  # 遇到错误立即退出

# ================================
# 🎨 颜色定义
# ================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color

# ================================
# 📊 全局变量
# ================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"  # plotumap 目录
CORE_DIR="$PROJECT_ROOT/core"
CONFIG_DIR="$PROJECT_ROOT/config"

# 使用配置脚本获取路径
get_config() {
    python3 "$SCRIPT_DIR/get_config.py" "$1"
}

# 从配置系统获取路径
BASE_DIR="$(get_config BASE)"
OUTPUT_DIR="$(get_config OUTPUT_DIR)"
IMPUTED_DIR="$(get_config IMPUTED_MATRICES)"
LOG_DIR="$OUTPUT_DIR/logs"

# 默认参数
MODE="test"
STAGES=""
MAX_CELLS=""
TEST_CELLS=3
USE_DISTANCE=true
SKIP_MERGE=false
SKIP_H5AD=false
SKIP_UMAP=false
GPU_TYPE="auto"
VERBOSE=true

# 运行统计
START_TIME=$(date +%s)
TOTAL_STAGES=0
SUCCESS_STAGES=0
FAILED_STAGES=0

# ================================
# 🔧 工具函数
# ================================

print_banner() {
    echo -e "${CYAN}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║                  🚀 HiRES Hi-C 处理器                        ║"
    echo "║                    一键启动脚本 v1.0                         ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
    mkdir -p "$LOG_DIR"  # 确保日志目录存在
    [[ $VERBOSE == true ]] && echo "$(date '+%Y-%m-%d %H:%M:%S') [INFO] $1" >> "$LOG_DIR/run.log"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
    mkdir -p "$LOG_DIR"  # 确保日志目录存在
    [[ $VERBOSE == true ]] && echo "$(date '+%Y-%m-%d %H:%M:%S') [WARN] $1" >> "$LOG_DIR/run.log"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
    mkdir -p "$LOG_DIR"  # 确保日志目录存在
    [[ $VERBOSE == true ]] && echo "$(date '+%Y-%m-%d %H:%M:%S') [ERROR] $1" >> "$LOG_DIR/run.log"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
    mkdir -p "$LOG_DIR"  # 确保日志目录存在
    [[ $VERBOSE == true ]] && echo "$(date '+%Y-%m-%d %H:%M:%S') [SUCCESS] $1" >> "$LOG_DIR/run.log"
}

show_progress() {
    local current=$1
    local total=$2
    local stage=$3
    local step=$4
    
    local percentage=$((current * 100 / total))
    local bar_length=30
    local filled_length=$((percentage * bar_length / 100))
    
    printf "\r${BLUE}[进度]${NC} [$stage] $step: ["
    printf "%*s" $filled_length | tr ' ' '█'
    printf "%*s" $((bar_length - filled_length)) | tr ' ' '░'
    printf "] %d%% (%d/%d)" $percentage $current $total
}

# ================================
# 🔍 环境检测
# ================================

detect_gpu() {
    if [[ $(uname -m) == "arm64" ]] && python3 -c "import torch; print(torch.backends.mps.is_available())" 2>/dev/null | grep -q "True"; then
        echo "mps"
    elif command -v nvidia-smi >/dev/null 2>&1; then
        echo "cuda"
    else
        echo "cpu"
    fi
}

check_environment() {
    log_info "🔍 检查运行环境..."
    
    # 检查Python环境
    if ! command -v python3 >/dev/null 2>&1; then
        log_error "未找到 Python3，请先安装 Python"
        exit 1
    fi
    
    # 检查必要目录
    local required_dirs=("$CORE_DIR" "$CONFIG_DIR" "$IMPUTED_DIR")
    for dir in "${required_dirs[@]}"; do
        if [[ ! -d "$dir" ]]; then
            log_error "必要目录不存在: $dir"
            exit 1
        fi
    done
    
    # 创建输出目录
    mkdir -p "$OUTPUT_DIR"/{merged_matrices,hic_h5ad_files,umap_plots,logs,reports}
    
    # 检测GPU类型
    if [[ $GPU_TYPE == "auto" ]]; then
        GPU_TYPE=$(detect_gpu)
    fi
    
    log_info "✅ 环境检查完成"
    log_info "📱 检测到GPU类型: $GPU_TYPE"
    log_info "📁 工作目录: $BASE_DIR"
    log_info "🐍 Python版本: $(python3 --version)"
}

# ================================
# 🔎 Stage 发现
# ================================

discover_stages() {
    log_info "🔎 自动发现可用的 stages..."
    
    local available_stages=()
    
    # 扫描 imputed_matrices_by_stage 目录
    if [[ -d "$IMPUTED_DIR" ]]; then
        while IFS= read -r -d '' dir; do
            local stage_name=$(basename "$dir")
            # 检查是否有细胞数据
            local cell_count=$(find "$dir" -maxdepth 1 -type d -name "*E*" | wc -l)
            if [[ $cell_count -gt 0 ]]; then
                available_stages+=("$stage_name")
                log_info "  📋 发现 stage: $stage_name ($cell_count 个细胞)"
            fi
        done < <(find "$IMPUTED_DIR" -maxdepth 1 -type d -name "E*" -print0)
    fi
    
    if [[ ${#available_stages[@]} -eq 0 ]]; then
        log_error "未发现任何可用的 stage 数据"
        exit 1
    fi
    
    # 根据模式选择要处理的stages
    case $MODE in
        "all")
            STAGES="${available_stages[*]}"
            ;;
        "test")
            # 测试模式选择前3个stages
            local test_stages=("${available_stages[@]:0:3}")
            STAGES="${test_stages[*]}"
            ;;
        "stage")
            if [[ -z "$STAGES" ]]; then
                log_error "stage 模式需要指定 --stages 参数"
                exit 1
            fi
            ;;
        *)
            log_error "未知模式: $MODE"
            exit 1
            ;;
    esac
    
    # 转换为数组
    IFS=' ' read -ra STAGE_ARRAY <<< "$STAGES"
    TOTAL_STAGES=${#STAGE_ARRAY[@]}
    
    log_info "🎯 将处理 $TOTAL_STAGES 个 stages: ${STAGES// /, }"
}

# ================================
# 🔧 处理函数
# ================================

merge_matrices() {
    local stage=$1
    log_info "🔗 开始合并 stage $stage 的细胞矩阵..."
    
    cd "$CORE_DIR"
    
    if python3 -c "
import sys
import os
sys.path.insert(0, '$CORE_DIR')

from matrix_merger_robust import RobustHiCMatrixMerger
from config import CONFIG

stage = '$stage'
stage_dir = '$IMPUTED_DIR/$stage'
output_dir = '$(get_config OUTPUT_DIR)/merged_matrices'

try:
    merger = RobustHiCMatrixMerger(CONFIG.CHROM_SIZES)
    merged_files = merger.merge_stage_cells_robust(stage_dir, stage, output_dir)
    success_count = len(merged_files)
    
    # 计算失败数量（估算）
    all_cells = [d for d in os.listdir(stage_dir) if os.path.isdir(os.path.join(stage_dir, d))]
    total_cells = len(all_cells)
    failed_count = total_cells - success_count
    
    print(f'合并完成: {success_count} 成功, {failed_count} 失败')
    sys.exit(0 if failed_count == 0 else 1)
except Exception as e:
    print(f'合并失败: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)
    "; then
        log_success "✅ Stage $stage 矩阵合并完成"
        return 0
    else
        log_error "❌ Stage $stage 矩阵合并失败"
        return 1
    fi
}

build_h5ad() {
    local stage=$1
    log_info "🧬 开始构建 stage $stage 的 h5ad 文件..."
    
    local stage_dir="$IMPUTED_DIR/$stage"
    # 使用配置系统的基础路径
    local metadata_file="$(get_config BASE)/data/hires/GSE223917_HiRES_emb_metadata.xlsx"
    local output_file="$(get_config H5AD_OUTPUT)/${stage}_processed.h5ad"
    
    # 检查metadata文件是否存在
    if [[ ! -f "$metadata_file" ]]; then
        log_error "Metadata文件不存在: $metadata_file"
        return 1
    fi
    
    # 构建参数
    local args=(
        "--stage-dir" "$stage_dir"
        "--obs-xlsx" "$metadata_file"
        "--output" "$output_file"
    )
    
    # 添加可选参数
    if [[ -n "$MAX_CELLS" && "$MODE" != "test" ]]; then
        args+=("--max-cells" "$MAX_CELLS")
    elif [[ "$MODE" == "test" ]]; then
        args+=("--max-cells" "$TEST_CELLS")
    fi
    
    if [[ $USE_DISTANCE == true ]]; then
        args+=("--distance-features")
    fi
    
    # 选择处理脚本
    local script=""
    case $GPU_TYPE in
        "mps")
            script="build_stage_h5ad_mps.py"
            ;;
        "cuda")
            script="build_stage_h5ad_cuda.py"
            ;;
        "cpu")
            script="build_stage_h5ad_simple.py"
            ;;
        *)
            log_error "未知GPU类型: $GPU_TYPE"
            return 1
            ;;
    esac
    
    cd "$CORE_DIR"
    
    if python3 "$script" "${args[@]}"; then
        log_success "✅ Stage $stage h5ad 构建完成: $output_file"
        return 0
    else
        log_error "❌ Stage $stage h5ad 构建失败"
        return 1
    fi
}

generate_umap() {
    local stage=$1
    log_info "📊 开始生成 stage $stage 的 UMAP 图..."
    
    local h5ad_file="$(get_config H5AD_OUTPUT)/${stage}_processed.h5ad"
    local output_dir="$OUTPUT_DIR/umap_plots"
    
    if [[ ! -f "$h5ad_file" ]]; then
        log_error "h5ad 文件不存在: $h5ad_file"
        return 1
    fi
    
    cd "$CORE_DIR"
    
    if python3 -c "
import sys
import os
sys.path.insert(0, '$CORE_DIR')

from hic_umap_visualizer import HiCUMAPVisualizer
import anndata as ad

try:
    # 加载数据
    adata = ad.read_h5ad('$h5ad_file')
    
    # 创建可视化器 - 使用相对路径
    color_mapping_file = os.path.join('$CONFIG_DIR', 'color_mapping.json')
    visualizer = HiCUMAPVisualizer(color_mapping_file)
    
    # 生成UMAP
    output_path = '$output_dir/${stage}_umap.png'
    visualizer.plot_umap(adata, title='$stage Hi-C UMAP', output_path=output_path)
    
    print(f'UMAP 生成完成: {output_path}')
    sys.exit(0)
except Exception as e:
    print(f'UMAP 生成失败: {e}')
    sys.exit(1)
    "; then
        log_success "✅ Stage $stage UMAP 生成完成"
        return 0
    else
        log_error "❌ Stage $stage UMAP 生成失败"
        return 1
    fi
}

# ================================
# 🚀 主处理流程
# ================================

process_stage() {
    local stage=$1
    local stage_num=$2
    
    log_info ""
    log_info "🎯 处理 Stage $stage ($stage_num/$TOTAL_STAGES)"
    log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    local stage_success=true
    
    # 步骤1: 合并矩阵
    if [[ $SKIP_MERGE == false ]]; then
        show_progress $((stage_num * 3 - 2)) $((TOTAL_STAGES * 3)) "$stage" "合并矩阵"
        if ! merge_matrices "$stage"; then
            stage_success=false
        fi
        echo  # 换行
    fi
    
    # 步骤2: 构建h5ad
    if [[ $SKIP_H5AD == false && $stage_success == true ]]; then
        show_progress $((stage_num * 3 - 1)) $((TOTAL_STAGES * 3)) "$stage" "构建h5ad"
        if ! build_h5ad "$stage"; then
            stage_success=false
        fi
        echo  # 换行
    fi
    
    # 步骤3: 生成UMAP
    if [[ $SKIP_UMAP == false && $stage_success == true ]]; then
        show_progress $((stage_num * 3)) $((TOTAL_STAGES * 3)) "$stage" "生成UMAP"
        if ! generate_umap "$stage"; then
            stage_success=false
        fi
        echo  # 换行
    fi
    
    if [[ $stage_success == true ]]; then
        ((SUCCESS_STAGES++))
        log_success "🎉 Stage $stage 处理完成"
    else
        ((FAILED_STAGES++))
        log_error "💥 Stage $stage 处理失败"
    fi
}

run_pipeline() {
    log_info ""
    log_info "🚀 开始执行 Hi-C 数据处理流水线"
    log_info "📊 处理模式: $MODE"
    log_info "⚡ GPU类型: $GPU_TYPE"
    log_info "🎯 距离特征: $([ $USE_DISTANCE == true ] && echo '启用' || echo '禁用')"
    
    # 处理每个stage
    local stage_num=1
    for stage in "${STAGE_ARRAY[@]}"; do
        process_stage "$stage" $stage_num
        ((stage_num++))
    done
}

# ================================
# 📊 结果汇总
# ================================

generate_summary() {
    local end_time=$(date +%s)
    local duration=$((end_time - START_TIME))
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))
    
    log_info ""
    log_info "📊 处理完成 - 结果汇总"
    log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log_info "⏱️  总耗时: ${hours}h ${minutes}m ${seconds}s"
    log_info "📋 总stages: $TOTAL_STAGES"
    log_success "✅ 成功: $SUCCESS_STAGES"
    [[ $FAILED_STAGES -gt 0 ]] && log_error "❌ 失败: $FAILED_STAGES"
    log_info ""
    log_info "📁 输出文件位置:"
    log_info "   🔗 合并矩阵: $OUTPUT_DIR/merged_matrices"
    log_info "   🧬 H5AD文件: $OUTPUT_DIR/hic_h5ad_files"
    log_info "   📊 UMAP图表: $OUTPUT_DIR/umap_plots"
    log_info "   📝 运行日志: $LOG_DIR/run.log"
    
    # 生成JSON报告
    cat > "$OUTPUT_DIR/reports/summary.json" << EOF
{
    "timestamp": "$(date -Iseconds)",
    "mode": "$MODE",
    "gpu_type": "$GPU_TYPE",
    "total_stages": $TOTAL_STAGES,
    "success_stages": $SUCCESS_STAGES,
    "failed_stages": $FAILED_STAGES,
    "duration_seconds": $duration,
    "stages_processed": [$(printf '"%s",' "${STAGE_ARRAY[@]}" | sed 's/,$//')],
    "use_distance_features": $USE_DISTANCE
}
EOF
    
    log_info "📄 详细报告: $OUTPUT_DIR/reports/summary.json"
}

# ================================
# 📖 帮助信息
# ================================

show_help() {
    echo -e "${WHITE}🚀 HiRES Hi-C 数据处理一键启动脚本${NC}"
    echo
    echo -e "${CYAN}用法:${NC}"
    echo "  $0 [选项]"
    echo
    echo -e "${CYAN}运行模式:${NC}"
    echo -e "  ${GREEN}--mode test${NC}         快速功能测试 (3个细胞 × 前3个stages)"
    echo -e "  ${GREEN}--mode all${NC}          处理所有发现的stages"
    echo -e "  ${GREEN}--mode stage${NC}        处理指定的stages (需配合 --stages)"
    echo
    echo -e "${CYAN}参数选项:${NC}"
    echo -e "  ${YELLOW}--stages STAGES${NC}     指定处理的stages (逗号分隔，如: E70,E80,EX05)"
    echo -e "  ${YELLOW}--max-cells N${NC}       每个stage最大处理细胞数 (默认: 全部)"
    echo -e "  ${YELLOW}--test-cells N${NC}      测试模式下的细胞数 (默认: 3)"
    echo -e "  ${YELLOW}--gpu-type TYPE${NC}     GPU类型 [auto|mps|cuda|cpu] (默认: auto)"
    echo
    echo -e "${CYAN}处理控制:${NC}"
    echo -e "  ${PURPLE}--use-distance${NC}      使用距离特征 (更快，推荐)"
    echo -e "  ${PURPLE}--skip-merge${NC}        跳过矩阵合并步骤"
    echo -e "  ${PURPLE}--skip-h5ad${NC}         跳过h5ad构建步骤"
    echo -e "  ${PURPLE}--skip-umap${NC}         跳过UMAP可视化步骤"
    echo
    echo -e "${CYAN}示例用法:${NC}"
    echo -e "  ${GREEN}$0 --mode test${NC}                           # 快速功能测试"
    echo -e "  ${GREEN}$0 --mode all --use-distance${NC}            # 处理所有stages（推荐）"
    echo -e "  ${GREEN}$0 --mode stage --stages E70,E80${NC}        # 处理指定stages"
    echo -e "  ${GREEN}$0 --mode test --gpu-type cpu${NC}           # 强制使用CPU测试"
    echo
    echo -e "${CYAN}说明:${NC}"
    echo "  🔍 脚本会自动发现 /output/imputed_matrices_by_stage 中的stages"
    echo "  ⚡ 自动检测硬件类型并选择最优处理方式 (MPS/CUDA/CPU)"
    echo "  📊 支持完整流水线: 矩阵合并 → h5ad构建 → UMAP可视化"
    echo "  🎯 推荐使用 --use-distance 选项以获得更好的性能"
}

# ================================
# 🔧 参数解析
# ================================

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --mode)
                MODE="$2"
                shift 2
                ;;
            --stages)
                STAGES="${2//,/ }"  # 将逗号替换为空格
                shift 2
                ;;
            --max-cells)
                MAX_CELLS="$2"
                shift 2
                ;;
            --test-cells)
                TEST_CELLS="$2"
                shift 2
                ;;
            --gpu-type)
                GPU_TYPE="$2"
                shift 2
                ;;
            --use-distance)
                USE_DISTANCE=true
                shift
                ;;
            --skip-merge)
                SKIP_MERGE=true
                shift
                ;;
            --skip-h5ad)
                SKIP_H5AD=true
                shift
                ;;
            --skip-umap)
                SKIP_UMAP=true
                shift
                ;;
            --help|-h)
                show_help
                exit 0
                ;;
            *)
                log_error "未知参数: $1"
                echo "使用 --help 查看帮助信息"
                exit 1
                ;;
        esac
    done
    
    # 参数验证
    if [[ ! "$MODE" =~ ^(test|all|stage)$ ]]; then
        log_error "无效的模式: $MODE (支持: test, all, stage)"
        exit 1
    fi
    
    if [[ ! "$GPU_TYPE" =~ ^(auto|mps|cuda|cpu)$ ]]; then
        log_error "无效的GPU类型: $GPU_TYPE (支持: auto, mps, cuda, cpu)"
        exit 1
    fi
}

# ================================
# 🚀 主函数
# ================================

main() {
    # 解析参数
    parse_arguments "$@"
    
    # 显示横幅
    print_banner
    
    # 检查环境
    check_environment
    
    # 发现stages
    discover_stages
    
    # 运行流水线
    run_pipeline
    
    # 生成汇总
    generate_summary
    
    log_info ""
    if [[ $FAILED_STAGES -eq 0 ]]; then
        log_success "🎉 所有任务完成！"
        exit 0
    else
        log_warn "⚠️  部分任务失败，请检查日志"
        exit 1
    fi
}

# 运行主函数
main "$@"
