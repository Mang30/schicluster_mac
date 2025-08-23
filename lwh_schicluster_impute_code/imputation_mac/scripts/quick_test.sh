#!/bin/bash
# M4 MacBook MPS 插补快速测试脚本
# 在 micromamba schicluster 环境中运行

echo "🚀 M4 MacBook MPS 插补快速测试"
echo "💻 硬件: M4 MacBook 24GB"
echo "🐍 环境: micromamba schicluster"
echo "=" * 50

# 激活 schicluster 环境 (macOS 兼容)
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

micromamba activate schicluster

# 检查 Python 和环境
echo "🔍 检查 Python 环境..."
python --version
echo "📍 Python 路径: $(which python)"
echo "📦 当前环境: $CONDA_DEFAULT_ENV"

# 运行环境测试
echo "🧪 运行环境测试..."
cd /Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac

# 确保在schicluster环境中运行
eval "$(micromamba shell hook --shell zsh)"
micromamba activate schicluster

# 验证环境激活
if [[ "$CONDA_DEFAULT_ENV" != "schicluster" ]]; then
    echo "❌ schicluster 环境激活失败，请手动激活后重试"
    echo "   命令: micromamba activate schicluster"
    exit 1
fi

python tests/test_environment.py

# 检查测试结果
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ 环境测试通过! 可以开始插补测试"
    echo ""
    echo "📋 下一步选择:"
    echo "1. 🧪 单细胞测试: python ../core/batch_mps_imputation_v2.py --stages E70 --max-cells 1"
    echo "2. 📊 进度监控: python monitor_progress.py --once"
    echo "3. 🚀 完整任务: bash start_parallel_imputation.sh"
    echo ""
    
    read -p "是否运行单细胞测试? (y/N): " test_single
    if [[ $test_single =~ ^[Yy]$ ]]; then
        echo "🧪 开始单细胞测试..."
        # 确保在schicluster环境中运行
        eval "$(micromamba shell hook --shell zsh)"
        micromamba activate schicluster
        
        # 验证环境
        if [[ "$CONDA_DEFAULT_ENV" != "schicluster" ]]; then
            echo "❌ 环境激活失败，尝试手动激活..."
            echo "   请运行: micromamba activate schicluster"
            exit 1
        fi
        
        python core/batch_mps_imputation_v2.py --stages E70 --max-cells 1
        
        if [ $? -eq 0 ]; then
            echo "🎉 单细胞测试成功! 系统已准备好进行完整插补"
        else
            echo "❌ 单细胞测试失败，请检查错误信息"
        fi
    fi
else
    echo "❌ 环境测试失败，请解决问题后重试"
    exit 1
fi