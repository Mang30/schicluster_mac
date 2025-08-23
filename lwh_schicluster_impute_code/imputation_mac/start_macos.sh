#!/bin/bash
# macOS 优化版启动脚本 - 解决 micromamba 环境问题

echo "🍎 macOS 优化版 M4 MacBook MPS 插补系统"
echo "=" * 50

# 检查 micromamba 是否可用
if ! command -v micromamba &> /dev/null; then
    echo "❌ micromamba 未找到，请先安装 micromamba"
    exit 1
fi

# 获取当前 shell
CURRENT_SHELL=$(basename "$SHELL")
echo "🐚 当前 shell: $CURRENT_SHELL"

# 初始化 micromamba（兼容 zsh 和 bash）
echo "🔧 初始化 micromamba..."
if [[ "$CURRENT_SHELL" == "zsh" ]]; then
    eval "$(micromamba shell hook --shell zsh)"
elif [[ "$CURRENT_SHELL" == "bash" ]]; then
    eval "$(micromamba shell hook --shell bash)"
else
    echo "⚠️  未知 shell，尝试 bash 模式"
    eval "$(micromamba shell hook --shell bash)"
fi

# 激活环境
echo "🔧 激活 schicluster 环境..."
micromamba activate schicluster

# 验证环境
if ! python -c "import sys; print('Python:', sys.executable)" 2>/dev/null; then
    echo "❌ Python 环境有问题"
    exit 1
fi

echo "✅ 环境验证成功"
echo ""

# 显示菜单
echo "📋 选择操作:"
echo "1. 🧪 快速测试 (单细胞验证)"
echo "2. 🔧 环境检查"
echo "3. 🚀 开始插补 (需要额外确认)"
echo "4. 📊 查看进度"
echo "5. ❌ 退出"

read -p "请选择 (1-5): " choice

case $choice in
    1)
        echo "🧪 开始快速测试..."
        cd "$(dirname "${BASH_SOURCE[0]}")"
        python core/batch_mps_imputation_v2.py --stages E70 --max-cells 1
        ;;
    2)
        echo "🔧 运行环境检查..."
        cd "$(dirname "${BASH_SOURCE[0]}")"
        python tests/test_environment.py
        ;;
    3)
        echo "🚀 准备开始插补..."
        echo "⚠️  这将需要很长时间，确认继续?"
        read -p "(y/N): " confirm
        if [[ $confirm =~ ^[Yy]$ ]]; then
            cd "$(dirname "${BASH_SOURCE[0]}")"
            # 这里可以调用你的插补脚本
            echo "🚀 启动插补..."
            python core/batch_mps_imputation_v2.py --stages E70 E75 E80 E85 E95 EX05 EX15
        else
            echo "❌ 已取消"
        fi
        ;;
    4)
        echo "📊 查看进度..."
        # 这里可以添加进度查看逻辑
        echo "功能开发中..."
        ;;
    5)
        echo "👋 再见!"
        exit 0
        ;;
    *)
        echo "❌ 无效选择"
        exit 1
        ;;
esac
