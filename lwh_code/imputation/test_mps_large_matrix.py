#!/usr/bin/env python3
"""
测试 MPS 是否能处理大矩阵
"""

import sys
import logging
import torch

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 添加路径
sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')
sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES')

def test_mps_large_matrix():
    """测试 MPS 处理大矩阵的能力"""
    print("🧪 测试 MPS 处理大矩阵")
    
    try:
        from impute_bin_format import impute_bin_matrix
        
        # 测试文件
        test_file = "/Volumes/SumSung500/CSU/0_HiRES/converted_matrices_by_stage_20k/E70/GasdE701001/GasdE701001_chr1.txt"
        output_file = "/tmp/test_mps_output.npz"
        chrom_size_file = "/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes.txt"
        
        print(f"📁 测试文件: {test_file}")
        print(f"📤 输出文件: {output_file}")
        
        # 强制使用 MPS
        print(f"🚀 开始 MPS 测试...")
        success = impute_bin_matrix(test_file, output_file, chrom_size_file, 'chr1', use_mps=True)
        
        if success:
            print("✅ MPS 测试成功！")
            import os
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file) / (1024*1024)  # MB
                print(f"📊 输出文件大小: {file_size:.1f} MB")
        else:
            print("❌ MPS 测试失败")
            
    except Exception as e:
        print(f"❌ 测试出错: {e}")
        import traceback
        traceback.print_exc()

def test_mps_memory():
    """测试 MPS 内存情况"""
    print("\n🧠 测试 MPS 内存")
    
    if not torch.backends.mps.is_available():
        print("❌ MPS 不可用")
        return
    
    device = torch.device("mps")
    print(f"✅ MPS 设备可用: {device}")
    
    # 测试不同大小的矩阵
    sizes = [1000, 5000, 8000, 10000]
    
    for size in sizes:
        try:
            print(f"\n测试 {size}x{size} 矩阵:")
            
            # 创建测试矩阵
            A = torch.randn(size, size, device=device, dtype=torch.float32)
            print(f"  ✅ 创建矩阵成功")
            
            # 简单运算测试
            B = torch.mm(A, A)
            print(f"  ✅ 矩阵乘法成功")
            
            # 清理内存
            del A, B
            if hasattr(torch.mps, 'empty_cache'):
                torch.mps.empty_cache()
            print(f"  ✅ 内存清理完成")
            
        except Exception as e:
            print(f"  ❌ 失败: {e}")
            break

if __name__ == "__main__":
    # 测试 MPS 内存
    test_mps_memory()
    
    # 测试实际插补
    test_mps_large_matrix()
