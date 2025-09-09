#!/bin/bash

# JSON to H5AD 批处理脚本
# 处理所有7个发育阶段的数据
# 版本: v1.0
# 创建日期: 2025-09-09

set -e  # 出错时停止执行

# 配置
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_PATH="/Users/wuhaoliu/mamba/envs/schicluster/bin/python"
CONVERTER_SCRIPT="$SCRIPT_DIR/convert_json_to_h5ad.py"
CONFIG_FILE="$SCRIPT_DIR/config.yaml"
METADATA_FILE="$SCRIPT_DIR/GSE223917_HiRES_emb_metadata.xlsx"

# 日志文件
LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOG_DIR"
BATCH_LOG="$LOG_DIR/batch_processing_$(date +%Y%m%d_%H%M%S).log"

# 所有stages
STAGES=("E70" "E75" "E80" "E85" "E95" "EX05" "EX15")

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日志函数
log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$BATCH_LOG"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$BATCH_LOG"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$BATCH_LOG"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$BATCH_LOG"
}

# 检查依赖
check_dependencies() {
    log_info "检查依赖..."
    
    # 检查Python路径
    if [[ ! -f "$PYTHON_PATH" ]]; then
        log_error "Python路径不存在: $PYTHON_PATH"
        exit 1
    fi
    
    # 检查转换脚本
    if [[ ! -f "$CONVERTER_SCRIPT" ]]; then
        log_error "转换脚本不存在: $CONVERTER_SCRIPT"
        exit 1
    fi
    
    # 检查配置文件
    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_error "配置文件不存在: $CONFIG_FILE"
        exit 1
    fi
    
    # 检查元数据文件
    if [[ ! -f "$METADATA_FILE" ]]; then
        log_warning "元数据文件不存在: $METADATA_FILE"
        log_info "将使用配置文件中指定的元数据文件路径"
    fi
    
    # 检查Python包
    log_info "检查Python环境..."
    "$PYTHON_PATH" -c "import pandas, numpy, scanpy, yaml, openpyxl; print('所有必需包已安装')" || {
        log_error "缺少必需的Python包"
        exit 1
    }
    
    log_success "依赖检查通过"
}

# 处理单个stage
process_stage() {
    local stage=$1
    local debug_mode=${2:-false}
    
    log_info "开始处理 Stage $stage"
    
    # 构建命令
    local cmd="$PYTHON_PATH $CONVERTER_SCRIPT --stage $stage --config $CONFIG_FILE"
    
    if [[ -f "$METADATA_FILE" ]]; then
        cmd="$cmd --metadata $METADATA_FILE"
    fi
    
    if [[ "$debug_mode" == "true" ]]; then
        cmd="$cmd --debug"
    fi
    
    # 执行命令
    local start_time=$(date +%s)
    
    if eval "$cmd"; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_success "Stage $stage 处理完成 (用时: ${duration}s)"
        return 0
    else
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        log_error "Stage $stage 处理失败 (用时: ${duration}s)"
        return 1
    fi
}

# 检查磁盘空间
check_disk_space() {
    local output_dir
    output_dir=$("$PYTHON_PATH" -c "
import yaml
with open('$CONFIG_FILE', 'r') as f:
    config = yaml.safe_load(f)
print(config['paths']['output_dir'])
" 2>/dev/null) || {
        log_warning "无法读取输出目录配置，跳过磁盘空间检查"
        return 0
    }
    
    local available_space
    available_space=$(df -h "$output_dir" 2>/dev/null | awk 'NR==2 {print $4}' | sed 's/G//')
    
    if [[ -n "$available_space" && $(echo "$available_space < 5" | bc -l 2>/dev/null) == "1" ]]; then
        log_warning "输出目录可用空间不足5GB: ${available_space}GB"
    else
        log_info "磁盘空间检查通过"
    fi
}

# 生成处理报告
generate_report() {
    local success_count=$1
    local total_count=$2
    local failed_stages=("${@:3}")
    
    log_info "生成处理报告..."
    
    local report_file="$LOG_DIR/batch_processing_report_$(date +%Y%m%d_%H%M%S).txt"
    
    cat > "$report_file" << EOF
JSON to H5AD 批处理报告
=====================

处理时间: $(date)
总Stage数: $total_count
成功处理: $success_count
失败处理: $((total_count - success_count))

EOF

    if [[ ${#failed_stages[@]} -gt 0 ]]; then
        echo "失败的Stages:" >> "$report_file"
        for stage in "${failed_stages[@]}"; do
            echo "  - $stage" >> "$report_file"
        done
        echo "" >> "$report_file"
    fi
    
    # 添加输出文件信息
    local output_dir
    output_dir=$("$PYTHON_PATH" -c "
import yaml
with open('$CONFIG_FILE', 'r') as f:
    config = yaml.safe_load(f)
print(config['paths']['output_dir'])
" 2>/dev/null) || output_dir="unknown"
    
    echo "输出目录: $output_dir" >> "$report_file"
    echo "生成的文件:" >> "$report_file"
    
    if [[ -d "$output_dir" ]]; then
        find "$output_dir" -name "*.h5ad" -exec ls -lh {} \; | awk '{print "  - " $9 " (" $5 ")"}' >> "$report_file" 2>/dev/null || true
    fi
    
    log_info "处理报告已保存: $report_file"
}

# 显示帮助信息
show_help() {
    cat << EOF
用法: $0 [选项] [stages...]

选项:
  -h, --help      显示帮助信息
  -d, --debug     调试模式 (处理少量文件)
  -s, --serial    串行处理 (默认)
  -c, --continue  出错时继续处理其他stages
  --stages        指定要处理的stages (默认处理所有)

示例:
  $0                    # 处理所有stages
  $0 --debug            # 调试模式处理所有stages
  $0 E75 E80            # 只处理E75和E80
  $0 --continue         # 出错时继续处理
  
支持的stages: ${STAGES[*]}
EOF
}

# 主函数
main() {
    local debug_mode=false
    local continue_on_error=false
    local selected_stages=()
    
    # 解析命令行参数
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -d|--debug)
                debug_mode=true
                shift
                ;;
            -c|--continue)
                continue_on_error=true
                shift
                ;;
            -s|--serial)
                # 串行模式是默认的，这里只是为了兼容
                shift
                ;;
            --stages)
                shift
                while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                    selected_stages+=("$1")
                    shift
                done
                ;;
            E70|E75|E80|E85|E95|EX05|EX15)
                selected_stages+=("$1")
                shift
                ;;
            *)
                log_error "未知参数: $1"
                show_help
                exit 1
                ;;
        esac
    done
    
    # 如果没有指定stages，使用所有stages
    if [[ ${#selected_stages[@]} -eq 0 ]]; then
        selected_stages=("${STAGES[@]}")
    fi
    
    log_info "开始批处理 JSON to H5AD 转换"
    log_info "处理模式: $([ "$debug_mode" == "true" ] && echo "调试模式" || echo "正常模式")"
    log_info "出错策略: $([ "$continue_on_error" == "true" ] && echo "继续处理" || echo "停止处理")"
    log_info "要处理的stages: ${selected_stages[*]}"
    
    # 检查依赖
    check_dependencies
    
    # 检查磁盘空间
    check_disk_space
    
    # 处理stages
    local success_count=0
    local failed_stages=()
    local total_start_time=$(date +%s)
    
    for stage in "${selected_stages[@]}"; do
        if process_stage "$stage" "$debug_mode"; then
            ((success_count++))
        else
            failed_stages+=("$stage")
            if [[ "$continue_on_error" != "true" ]]; then
                log_error "Stage $stage 处理失败，停止批处理"
                break
            fi
        fi
        
        # 处理间隔 (避免系统负载过高)
        if [[ ${#selected_stages[@]} -gt 1 ]]; then
            log_info "等待5秒后处理下一个stage..."
            sleep 5
        fi
    done
    
    local total_end_time=$(date +%s)
    local total_duration=$((total_end_time - total_start_time))
    
    # 生成报告
    generate_report "$success_count" "${#selected_stages[@]}" "${failed_stages[@]}"
    
    # 输出最终结果
    log_info "批处理完成，总用时: ${total_duration}s"
    log_success "成功处理: $success_count/${#selected_stages[@]} 个stages"
    
    if [[ ${#failed_stages[@]} -gt 0 ]]; then
        log_error "失败的stages: ${failed_stages[*]}"
        exit 1
    else
        log_success "所有stages处理成功!"
        exit 0
    fi
}

# 捕获Ctrl+C信号
trap 'log_warning "用户中断批处理"; exit 130' INT

# 运行主函数
main "$@"