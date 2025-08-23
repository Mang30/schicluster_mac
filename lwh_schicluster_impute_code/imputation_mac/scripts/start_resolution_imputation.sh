#!/bin/bash
# M4 MacBook 24GB 一键分辨率调整并行插补启动脚本
# 支持动态分辨率设置

echo "🚀 M4 MacBook 24GB MPS 加速并行插补系统 (分辨率可调)"
echo "=" * 80
echo "💻 硬件配置: M4 MacBook 24GB 统一内存"
echo "🔧 加速技术: MPS (Metal Performance Shaders)"
echo "📊 预计加速: 10-15x vs CPU"
echo "🎯 分辨率支持: 5k, 10k, 20k, 40k, 100k, 200k, 500k, 1M"
echo "=" * 80

# 激活环境 (macOS 兼容)
echo "🔧 激活 schicluster 环境..."

# 检测 shell 类型并相应初始化 micromamba
if [[ "$SHELL" == *"zsh"* ]] || [[ -n "$ZSH_VERSION" ]]; then
    echo "📱 检测到 zsh shell (macOS 默认)"
    eval "$(micromamba shell hook --shell zsh)"
elif [[ "$SHELL" == *"bash"* ]] || [[ -n "$BASH_VERSION" ]]; then
    echo "🐚 检测到 bash shell"
    eval "$(micromamba shell hook --shell bash)"
else
    echo "⚠️  未知 shell，尝试通用初始化"
    eval "$(micromamba shell hook --shell bash)"
fi

# 激活环境
micromamba activate schicluster

# 检查环境
echo "🔍 检查环境依赖..."
python -c "import torch; print('✅ PyTorch 可用:', torch.__version__)" 2>/dev/null || echo "❌ PyTorch 未安装"
python -c "import torch; print('✅ MPS 可用:', torch.backends.mps.is_available())" 2>/dev/null || echo "❌ MPS 不可用"

# 分辨率选择菜单
echo ""
echo "🎯 请选择分辨率 (分辨率越高，计算时间越长，内存需求越大):"
echo "1. 5k   (5,000 bp)    - 超高分辨率 🔬"
echo "2. 10k  (10,000 bp)   - 高分辨率 📊"
echo "3. 20k  (20,000 bp)   - 中高分辨率 ⚡"
echo "4. 40k  (40,000 bp)   - 中等分辨率 🎯"
echo "5. 100k (100,000 bp)  - 标准分辨率 ✅ (推荐)"
echo "6. 200k (200,000 bp)  - 低分辨率 🚀"
echo "7. 500k (500,000 bp)  - 超低分辨率 ⚡⚡"
echo "8. 1M   (1,000,000 bp) - 极低分辨率 🚄"
echo "9. 🔧 自定义分辨率"

read -p "请选择分辨率 (1-9): " res_choice

case $res_choice in
    1)
        RESOLUTION=5000
        echo "🔬 选择: 5k 超高分辨率 (计算密集型)"
        ;;
    2)
        RESOLUTION=10000
        echo "📊 选择: 10k 高分辨率"
        ;;
    3)
        RESOLUTION=20000
        echo "⚡ 选择: 20k 中高分辨率"
        ;;
    4)
        RESOLUTION=40000
        echo "🎯 选择: 40k 中等分辨率"
        ;;
    5)
        RESOLUTION=100000
        echo "✅ 选择: 100k 标准分辨率 (推荐)"
        ;;
    6)
        RESOLUTION=200000
        echo "🚀 选择: 200k 低分辨率"
        ;;
    7)
        RESOLUTION=500000
        echo "⚡⚡ 选择: 500k 超低分辨率"
        ;;
    8)
        RESOLUTION=1000000
        echo "🚄 选择: 1M 极低分辨率"
        ;;
    9)
        read -p "请输入自定义分辨率 (例如: 50000): " RESOLUTION
        echo "🔧 选择: 自定义分辨率 ${RESOLUTION}bp"
        ;;
    *)
        echo "❌ 无效选择，使用默认 100k 分辨率"
        RESOLUTION=100000
        ;;
esac

# 验证分辨率是否为有效数字
if ! [[ "$RESOLUTION" =~ ^[0-9]+$ ]]; then
    echo "❌ 分辨率必须是数字，使用默认 100k"
    RESOLUTION=100000
fi

echo "📏 确认分辨率: ${RESOLUTION} bp"

# 显示执行模式菜单
echo ""
echo "📋 请选择执行模式:"
echo "1. 🧪 测试模式 (每个阶段处理2个细胞)"
echo "2. 🚀 完整并行执行 (所有6个任务)"
echo "3. 🎯 自定义任务选择"
echo "4. 📊 进度监控 (查看当前状态)"
echo "5. 📄 生成汇总报告"
echo "6. 🔧 单个任务执行"
echo "7. ❌ 退出"

read -p "请输入选择 (1-7): " choice

case $choice in
    1)
        echo "🧪 启动测试模式 (分辨率: ${RESOLUTION}bp)..."
        cd "$(dirname "${BASH_SOURCE[0]}")/../core"
        python batch_mps_imputation_v2.py \
            --stages E70 E75 \
            --max-cells 2 \
            --batch-size 2 \
            --resolution $RESOLUTION
        ;;
    2)
        echo "🚀 启动完整并行执行 (分辨率: ${RESOLUTION}bp)..."
        # 生成动态任务脚本
        generate_dynamic_tasks $RESOLUTION
        # 启动并行任务
        start_parallel_tasks
        ;;
    3)
        echo "🎯 自定义任务选择..."
        custom_task_selection $RESOLUTION
        ;;
    4)
        echo "📊 查看进度监控..."
        monitor_progress
        ;;
    5)
        echo "📄 生成汇总报告..."
        generate_summary_report
        ;;
    6)
        echo "🔧 单个任务执行..."
        single_task_execution $RESOLUTION
        ;;
    7)
        echo "❌ 退出"
        exit 0
        ;;
    *)
        echo "❌ 无效选择"
        exit 1
        ;;
esac

# 生成动态任务脚本函数
generate_dynamic_tasks() {
    local resolution=$1
    echo "🔧 生成分辨率 ${resolution}bp 的动态任务脚本..."
    
    # 创建临时任务目录
    TEMP_TASK_DIR="/tmp/hires_tasks_${resolution}"
    mkdir -p "$TEMP_TASK_DIR"
    
    # 生成6个任务脚本
    for i in {1..6}; do
        case $i in
            1) stages="E70";;
            2) stages="E75";;
            3) stages="E80";;
            4) stages="E85";;
            5) stages="E95";;
            6) stages="EX05 EX15";;
        esac
        
        cat > "${TEMP_TASK_DIR}/task${i}_imputation_${resolution}.sh" << EOF
#!/bin/bash
# 动态生成的任务${i} - 分辨率: ${resolution}bp
# M4 MacBook 24GB 优化

TASK_ID="task${i}"
STAGES="${stages}"
BATCH_SIZE=10
RESOLUTION=${resolution}

echo "🚀 开始任务 \$TASK_ID (分辨率: \${RESOLUTION}bp)"
echo "📦 处理阶段: \$STAGES"
echo "🔧 批次大小: \$BATCH_SIZE"

cd "$(dirname "\${BASH_SOURCE[0]}")/../core"

echo "🔧 执行命令: python batch_mps_imputation_v2.py --stages \$STAGES --batch-size \$BATCH_SIZE --resolution \$RESOLUTION"

python batch_mps_imputation_v2.py \\
    --stages \$STAGES \\
    --batch-size \$BATCH_SIZE \\
    --resolution \$RESOLUTION \\
    > "../logs/task${i}_\${RESOLUTION}bp_\$(date +%Y%m%d_%H%M%S).log" 2>&1

echo "✅ 任务 \$TASK_ID 完成 (分辨率: \${RESOLUTION}bp)"
EOF
        chmod +x "${TEMP_TASK_DIR}/task${i}_imputation_${resolution}.sh"
    done
    
    echo "✅ 动态任务脚本生成完成: $TEMP_TASK_DIR"
}

# 启动并行任务函数
start_parallel_tasks() {
    echo "🚀 启动并行任务..."
    
    # 创建日志目录
    mkdir -p "../logs"
    
    # 并行启动任务
    for i in {1..6}; do
        echo "🔥 启动任务 $i (分辨率: ${RESOLUTION}bp)..."
        "${TEMP_TASK_DIR}/task${i}_imputation_${RESOLUTION}.sh" &
        sleep 2  # 错开启动时间
    done
    
    echo "⏳ 所有任务已启动，等待完成..."
    wait
    echo "🎉 所有并行任务完成！"
}

# 自定义任务选择函数
custom_task_selection() {
    local resolution=$1
    echo "🎯 自定义任务选择 (分辨率: ${resolution}bp)"
    echo "请选择要执行的阶段 (可多选，用空格分隔):"
    echo "1. E70   2. E75   3. E80   4. E85   5. E95   6. EX05   7. EX15"
    
    read -p "请输入选择 (例如: 1 2 5): " stage_choices
    
    selected_stages=""
    for choice in $stage_choices; do
        case $choice in
            1) selected_stages="$selected_stages E70";;
            2) selected_stages="$selected_stages E75";;
            3) selected_stages="$selected_stages E80";;
            4) selected_stages="$selected_stages E85";;
            5) selected_stages="$selected_stages E95";;
            6) selected_stages="$selected_stages EX05";;
            7) selected_stages="$selected_stages EX15";;
        esac
    done
    
    if [ -z "$selected_stages" ]; then
        echo "❌ 未选择任何阶段"
        exit 1
    fi
    
    echo "✅ 选择的阶段: $selected_stages"
    read -p "批次大小 (默认10): " batch_size
    batch_size=${batch_size:-10}
    
    cd "$(dirname "${BASH_SOURCE[0]}")/../core"
    python batch_mps_imputation_v2.py \
        --stages $selected_stages \
        --batch-size $batch_size \
        --resolution $resolution
}

# 单个任务执行函数
single_task_execution() {
    local resolution=$1
    echo "🔧 单个任务执行 (分辨率: ${resolution}bp)"
    
    echo "请选择要执行的单个阶段:"
    echo "1. E70   2. E75   3. E80   4. E85   5. E95   6. EX05   7. EX15"
    read -p "请选择 (1-7): " stage_choice
    
    case $stage_choice in
        1) stage="E70";;
        2) stage="E75";;
        3) stage="E80";;
        4) stage="E85";;
        5) stage="E95";;
        6) stage="EX05";;
        7) stage="EX15";;
        *) echo "❌ 无效选择"; exit 1;;
    esac
    
    read -p "批次大小 (默认10): " batch_size
    batch_size=${batch_size:-10}
    
    read -p "最大细胞数 (默认无限制): " max_cells
    max_cells_arg=""
    if [ ! -z "$max_cells" ]; then
        max_cells_arg="--max-cells $max_cells"
    fi
    
    cd "$(dirname "${BASH_SOURCE[0]}")/../core"
    python batch_mps_imputation_v2.py \
        --stages $stage \
        --batch-size $batch_size \
        --resolution $resolution \
        $max_cells_arg
}

# 进度监控函数
monitor_progress() {
    echo "📊 系统状态监控"
    echo "=" * 50
    
    # GPU 使用情况
    echo "🔥 GPU 状态:"
    python -c "
import torch
if torch.backends.mps.is_available():
    print(f'✅ MPS 可用')
    print(f'💾 GPU 内存: 24GB 统一内存')
else:
    print('❌ MPS 不可用')
" 2>/dev/null
    
    # 进程状态
    echo ""
    echo "🔍 Python 进程:"
    ps aux | grep python | grep -v grep | head -5
    
    # 日志状态
    echo ""
    echo "📄 最新日志:"
    ls -lt ../logs/*.log 2>/dev/null | head -3 || echo "暂无日志文件"
    
    # 输出目录状态
    echo ""
    echo "📁 输出目录状态:"
    OUTPUT_DIR="/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage"
    if [ -d "$OUTPUT_DIR" ]; then
        echo "输出目录: $OUTPUT_DIR"
        find "$OUTPUT_DIR" -name "*.npz" | wc -l | xargs echo "已完成矩阵数:"
        du -sh "$OUTPUT_DIR" 2>/dev/null | cut -f1 | xargs echo "占用空间:"
    else
        echo "输出目录不存在"
    fi
}

# 汇总报告函数
generate_summary_report() {
    echo "📄 生成汇总报告..."
    
    REPORT_FILE="../logs/summary_report_$(date +%Y%m%d_%H%M%S).txt"
    
    {
        echo "🎯 Hi-C 插补汇总报告"
        echo "生成时间: $(date)"
        echo "=" * 60
        
        echo ""
        echo "📊 系统信息:"
        echo "硬件: M4 MacBook 24GB"
        echo "加速: MPS (Metal Performance Shaders)"
        uname -a
        
        echo ""
        echo "📁 输出统计:"
        OUTPUT_DIR="/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage"
        if [ -d "$OUTPUT_DIR" ]; then
            find "$OUTPUT_DIR" -name "*.npz" | wc -l | xargs echo "已完成矩阵数:"
            du -sh "$OUTPUT_DIR" 2>/dev/null | cut -f1 | xargs echo "占用空间:"
            
            echo ""
            echo "📦 各阶段进度:"
            for stage in E70 E75 E80 E85 E95 EX05 EX15; do
                if [ -d "$OUTPUT_DIR/$stage" ]; then
                    count=$(find "$OUTPUT_DIR/$stage" -name "*.npz" 2>/dev/null | wc -l)
                    echo "$stage: $count 个矩阵"
                fi
            done
        fi
        
        echo ""
        echo "📄 最新日志文件:"
        ls -lt ../logs/*.log 2>/dev/null | head -5 || echo "暂无日志文件"
        
    } > "$REPORT_FILE"
    
    echo "✅ 汇总报告生成完成: $REPORT_FILE"
    echo "📖 查看报告: cat $REPORT_FILE"
}

echo "🎉 脚本执行完成！"
