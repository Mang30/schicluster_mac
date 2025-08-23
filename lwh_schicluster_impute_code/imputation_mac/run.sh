#!/bin/bash
# M4 MacBook MPS 插补系统主启动脚本

echo "🚀 M4 MacBook MPS 加速单细胞 Hi-C 插补系统"
echo "=" * 70
echo "💻 硬件: M4 MacBook 24GB"
echo "🔧 MPS 加速: 10-15x vs CPU"
echo "📊 预计数据: ~7,476个细胞 × 20染色体 = ~149,520个矩阵"
echo "⏱️  预估时间: 12-20小时 (6任务并行)"
echo "=" * 70

# 检查当前目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "📍 当前目录: $SCRIPT_DIR"

# 检查目录结构
if [ ! -d "core" ] || [ ! -d "scripts" ] || [ ! -d "tests" ]; then
    echo "❌ 目录结构不完整，请确认在正确的项目根目录"
    exit 1
fi

echo ""
echo "📋 选择操作模式:"
echo "1. 🧪 快速测试 (环境验证 + 单细胞测试)"
echo "2. 🚀 完整并行执行 (6任务并行插补)"
echo "3. 📊 进度监控 (查看当前状态)"
echo "4. 🔧 环境测试 (检查MPS和依赖)"
echo "5. 📄 生成报告 (汇总执行结果)"
echo "6. ❌ 退出"

read -p "请输入选择 (1-6): " choice

case $choice in
    1)
        echo "🧪 启动快速测试..."
        bash scripts/quick_test.sh
        ;;
    2)
        echo "🚀 启动完整并行执行..."
        echo ""
        echo "⚠️  注意事项:"
        echo "   • 这将启动6个并行Python进程"
        echo "   • 预计需要12-20小时完成"
        echo "   • 将占用大量CPU和内存资源"
        echo "   • 建议在充电状态下运行"
        echo ""
        read -p "确认开始完整执行? (y/N): " confirm
        if [[ $confirm =~ ^[Yy]$ ]]; then
            bash scripts/start_parallel_imputation.sh
        else
            echo "❌ 操作已取消"
        fi
        ;;
    3)
        echo "📊 启动进度监控..."
        python core/monitor_progress.py
        ;;
    4)
        echo "🔧 运行环境测试..."
        python tests/test_environment.py
        ;;
    5)
        echo "📄 生成汇总报告..."
        python core/monitor_progress.py --summary
        ;;
    6)
        echo "👋 再见!"
        exit 0
        ;;
    *)
        echo "❌ 无效选择"
        exit 1
        ;;
esac

echo ""
echo "📋 后续可用命令:"
echo "   • 监控进度: python core/monitor_progress.py"
echo "   • 查看日志: tail -f logs/*.log"
echo "   • 系统监控: htop 或 Activity Monitor"
echo "   • 重新运行: bash run.sh"