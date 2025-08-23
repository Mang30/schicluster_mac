#!/bin/bash
# M4 MacBook 24GB 并行插补启动脚本
# 简化的启动界面和选项

echo "🚀 M4 MacBook 24GB MPS 加速并行插补系统"
echo "=" * 60
echo "💻 硬件配置: M4 MacBook 24GB 统一内存"
echo "🔧 加速技术: MPS (Metal Performance Shaders)"
echo "📊 预计加速: 10-15x vs CPU"
echo "=" * 60

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

# 生成动态任务脚本函数
generate_dynamic_tasks() {
    local resolution=$1
    echo "🔧 生成分辨率 ${resolution}bp 的动态任务脚本..."
    
    # 创建临时任务目录
    TEMP_TASK_DIR="/tmp/hires_tasks_${resolution}"
    mkdir -p "$TEMP_TASK_DIR"
    
    # 生成6个任务脚本 - 新的任务分配方案
    for i in {1..6}; do
        case $i in
            1) stages="E70 E80";;     # 任务1：E70 和 E80
            2) stages="E75";;         # 任务2：E75 独占
            3) stages="E85";;         # 任务3：E85 独占
            4) stages="E95";;         # 任务4：E95 独占
            5) stages="EX05";;        # 任务5：EX05 独占
            6) stages="EX15";;        # 任务6：EX15 独占
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

cd "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac/core"

echo "🔧 执行命令: python batch_mps_imputation_v2.py --stages \$STAGES --batch-size \$BATCH_SIZE --resolution \$RESOLUTION"

python batch_mps_imputation_v2.py \\
    --stages \$STAGES \\
    --batch-size \$BATCH_SIZE \\
    --resolution \$RESOLUTION \\
    > "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac/logs/task${i}_\${RESOLUTION}bp_\$(date +%Y%m%d_%H%M%S).log" 2>&1

echo "✅ 任务 \$TASK_ID 完成 (分辨率: \${RESOLUTION}bp)"
EOF
        chmod +x "${TEMP_TASK_DIR}/task${i}_imputation_${resolution}.sh"
    done
    
    echo "✅ 动态任务脚本生成完成: $TEMP_TASK_DIR"
}

# 启动并行任务函数
start_parallel_tasks() {
    local resolution=$1
    echo "🚀 启动并行任务 (分辨率: ${resolution}bp)..."
    
    # 创建日志目录（使用绝对路径）
    mkdir -p "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac/logs"
    
    # 并行启动任务
    for i in {1..6}; do
        echo "🔥 启动任务 $i (分辨率: ${resolution}bp)..."
        "${TEMP_TASK_DIR}/task${i}_imputation_${resolution}.sh" &
        sleep 2  # 错开启动时间
    done
    
    echo "⏳ 所有任务已启动，等待完成..."
    wait
    echo "🎉 所有并行任务完成！"
}

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

read -p "请选择分辨率 (1-9, 默认5): " res_choice
res_choice=${res_choice:-5}

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

# 显示菜单
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
        cd "$(dirname "${BASH_SOURCE[0]}")/.."
        python core/batch_mps_imputation_v2.py --stages E70 --max-cells 2 --resolution $RESOLUTION
        ;;
    2)
        echo "🚀 启动完整并行执行 (分辨率: ${RESOLUTION}bp)..."
        echo ""
        echo "📋 任务分配方案:"
        echo "   • 任务1: E70 + E80 (557 + 559 = 1,116 细胞)"
        echo "   • 任务2: E75 (1,870 细胞)"
        echo "   • 任务3: E85 (753 细胞)" 
        echo "   • 任务4: E95 (1,108 细胞)"
        echo "   • 任务5: EX05 (916 细胞)"
        echo "   • 任务6: EX15 (1,122 细胞)"
        echo ""
        echo "⚠️  注意: 这将启动6个并行任务，预计需要12-20小时"
        read -p "确认继续? (y/N): " confirm
        if [[ $confirm =~ ^[Yy]$ ]]; then
            # 生成动态任务脚本并启动
            generate_dynamic_tasks $RESOLUTION
            start_parallel_tasks $RESOLUTION
        else
            echo "❌ 操作已取消"
        fi
        ;;
    3)
        echo "🎯 自定义任务选择 (分辨率: ${RESOLUTION}bp):"
        echo "选择要处理的 stages:"
        echo "1. E70 (557 细胞)"
        echo "2. E75 (1,870 细胞)"  
        echo "3. E80 (559 细胞)"
        echo "4. E85 (753 细胞)"
        echo "5. E95 (1,108 细胞)"
        echo "6. EX05 (916 细胞)"
        echo "7. EX15 (1,122 细胞)"
        read -p "输入要处理的 stages (例: E70 E75 E80): " stages
        if [[ -n "$stages" ]]; then
            cd "$(dirname "${BASH_SOURCE[0]}")/.."
            python core/batch_mps_imputation_v2.py --stages $stages --resolution $RESOLUTION
        else
            echo "❌ 未输入有效的 stages"
        fi
        ;;
    4)
        echo "📊 启动进度监控..."
        echo "功能开发中... 可以查看日志文件了解进度"
        ls -la logs/ 2>/dev/null || echo "暂无日志文件"
        ;;
    5)
        echo "📄 生成汇总报告..."
        echo "功能开发中... 可以手动检查输出目录"
        echo "输出目录: /Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage/"
        ;;
    6)
        echo "🔧 单个 stage 执行:"
        echo "1. E70 (557 细胞)"
        echo "2. E75 (1,870 细胞)"
        echo "3. E80 (559 细胞)" 
        echo "4. E85 (753 细胞)"
        echo "5. E95 (1,108 细胞)"
        echo "6. EX05 (916 细胞)"
        echo "7. EX15 (1,122 细胞)"
        read -p "选择 stage (1-7): " stage_num
        case $stage_num in
            1) stage="E70" ;;
            2) stage="E75" ;;
            3) stage="E80" ;;
            4) stage="E85" ;;
            5) stage="E95" ;;
            6) stage="EX05" ;;
            7) stage="EX15" ;;
            *) echo "❌ 无效选择"; exit 1 ;;
        esac
        echo "🚀 启动 $stage (分辨率: ${RESOLUTION}bp)..."
        cd "$(dirname "${BASH_SOURCE[0]}")/.."
        python core/batch_mps_imputation_v2.py --stages $stage --resolution $RESOLUTION
        ;;
    7)
        echo "👋 再见!"
        exit 0
        ;;
    *)
        echo "❌ 无效选择，请重试"
        exit 1
        ;;
esac

# 显示后续操作提示
echo ""
echo "📋 后续操作提示:"
echo "• 查看输出: ls /Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage/"
echo "• 查看日志: tail -f ../logs/*.log"  
echo "• 系统监控: htop 或 Activity Monitor"
echo "• 停止任务: Ctrl+C"