#!/usr/bin/env python3
"""
快速验证脚本 - 验证 Hi-C NPZ 到 AnnData 转换工具的基本功能

运行这个脚本来快速检查工具是否能正常工作。
"""

import sys
from pathlib import Path
import logging

# 添加当前目录到 Python 路径
sys.path.insert(0, str(Path(__file__).parent))

try:
    from hic_converter import (
        HiCNpzLoader, 
        UpperTriangleExtractor, 
        ChromosomeManager,
        FeatureGenerator, 
        AnnDataBuilder
    )
    print("✅ 成功导入所有核心模块")
except ImportError as e:
    print(f"❌ 导入失败: {e}")
    sys.exit(1)

def check_dependencies():
    """检查依赖包"""
    required_packages = [
        'numpy', 'scipy', 'pandas', 'anndata'
    ]
    
    missing = []
    for package in required_packages:
        try:
            __import__(package)
            print(f"✅ {package}")
        except ImportError:
            print(f"❌ {package} - 未安装")
            missing.append(package)
    
    return missing

def check_input_data():
    """检查输入数据结构"""
    base_dir = Path("../output/imputed_matrices_by_stage")
    
    if not base_dir.exists():
        print(f"❌ 输入目录不存在: {base_dir}")
        return False
    
    stages = [d for d in base_dir.iterdir() if d.is_dir() and d.name.startswith(('E', 'EX'))]
    if not stages:
        print(f"❌ 在 {base_dir} 中没有找到发育阶段目录")
        return False
    
    print(f"✅ 找到 {len(stages)} 个发育阶段: {[s.name for s in stages[:3]]}...")
    
    # 检查第一个阶段的细胞数据
    first_stage = stages[0]
    cells = [d for d in first_stage.iterdir() if d.is_dir()]
    
    if not cells:
        print(f"❌ 在阶段 {first_stage.name} 中没有找到细胞目录")
        return False
    
    print(f"✅ 在阶段 {first_stage.name} 中找到 {len(cells)} 个细胞")
    
    # 检查染色体文件
    first_cell = cells[0]
    chr_files = list(first_cell.glob("*_chr*_*.npz"))
    
    if not chr_files:
        print(f"❌ 在细胞 {first_cell.name} 中没有找到染色体文件")
        return False
    
    print(f"✅ 在细胞 {first_cell.name} 中找到 {len(chr_files)} 个染色体文件")
    
    return True

def main():
    """主验证函数"""
    print("🔍 Hi-C NPZ 到 AnnData 转换工具验证")
    print("=" * 50)
    
    # 检查依赖
    print("\n📦 检查Python依赖包:")
    missing = check_dependencies()
    
    if missing:
        print(f"\n❌ 缺少依赖包: {missing}")
        print("请安装缺少的包:")
        print(f"pip install {' '.join(missing)}")
        return 1
    
    # 检查输入数据
    print("\n📁 检查输入数据结构:")
    if not check_input_data():
        print("\n❌ 输入数据检查失败")
        print("请确保数据位于 ../output/imputed_matrices_by_stage/ 目录")
        return 1
    
    # 检查组件初始化
    print("\n🔧 检查工具组件:")
    try:
        loader = HiCNpzLoader()
        extractor = UpperTriangleExtractor(include_diagonal=True)
        chr_manager = ChromosomeManager()
        feature_generator = FeatureGenerator(loader, extractor, chr_manager)
        anndata_builder = AnnDataBuilder()
        print("✅ 所有组件初始化成功")
    except Exception as e:
        print(f"❌ 组件初始化失败: {e}")
        return 1
    
    print("\n🎉 所有检查通过！")
    print("💡 现在可以运行:")
    print("   python test_conversion.py      # 运行完整测试")
    print("   python convert_npz_to_h5ad.py  # 开始转换")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
