#!/bin/bash
# 检查schicluster环境中的依赖

echo "🔍 检查 schicluster 环境依赖"
echo "=" * 50

# 激活环境
eval "$(micromamba shell hook --shell bash)"
micromamba activate schicluster

# 检查Python和环境
echo "📍 Python路径: $(which python)"
echo "📦 当前环境: $CONDA_DEFAULT_ENV"

# 检查关键包
echo ""
echo "🔧 检查关键依赖包:"

packages=("torch" "cooler" "scipy" "numpy" "pandas" "psutil")

for package in "${packages[@]}"; do
    if python -c "import $package" 2>/dev/null; then
        version=$(python -c "import $package; print($package.__version__)" 2>/dev/null || echo "unknown")
        echo "✅ $package: $version"
    else
        echo "❌ $package: 未安装"
    fi
done

# 检查MPS
echo ""
echo "🎯 检查MPS支持:"
python -c "
try:
    import torch
    print('✅ PyTorch 已安装')
    if torch.backends.mps.is_available():
        print('✅ MPS 可用')
        device = torch.device('mps')
        test_tensor = torch.randn(10, 10, device=device)
        result = torch.mm(test_tensor, test_tensor.t())
        print('✅ MPS 设备测试成功')
    else:
        print('⚠️  MPS 不可用')
except ImportError:
    print('❌ PyTorch 未安装')
except Exception as e:
    print(f'❌ MPS 测试失败: {e}')
"