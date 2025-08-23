#!/usr/bin/env python3
"""
测试 MPS GPU 加速插补功能的简单脚本
"""

import sys
import logging
import numpy as np
from scipy.sparse import csr_matrix

# 设置日志级别
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_mps_availability():
    """测试 MPS 是否可用"""
    print("=== MPS 可用性测试 ===")
    
    try:
        import torch
        print(f"✅ PyTorch 版本: {torch.__version__}")
        
        if torch.backends.mps.is_available():
            print("✅ MPS 加速可用")
            return True
        else:
            print("❌ MPS 加速不可用")
            return False
    except ImportError:
        print("❌ PyTorch 未安装")
        return False

def test_impute_functions():
    """测试插补函数"""
    print("\n=== 插补函数测试 ===")
    
    try:
        # 导入修改后的函数
        sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')
        from schicluster.impute.impute_chromosome import (
            random_walk_cpu, 
            random_walk_auto,
            MPS_AVAILABLE
        )
        
        print("✅ 成功导入插补函数")
        
        # 创建测试矩阵
        size = 100
        np.random.seed(42)
        data = np.random.random(size * 10)
        row = np.random.randint(0, size, size * 10)
        col = np.random.randint(0, size, size * 10)
        P = csr_matrix((data, (row, col)), shape=(size, size))
        
        # 归一化为概率矩阵
        P = P + P.T  # 对称化
        row_sums = np.array(P.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1  # 避免除零
        P = P.multiply(1/row_sums[:, np.newaxis])
        
        print(f"✅ 创建测试矩阵: {size}x{size}, 稀疏度: {P.nnz/(size*size):.4f}")
        
        # 测试 CPU 版本
        import time
        start_time = time.time()
        result_cpu = random_walk_cpu(P, rp=0.5, tol=0.01)
        cpu_time = time.time() - start_time
        print(f"✅ CPU 版本完成，耗时: {cpu_time:.3f} 秒")
        
        # 测试自动选择版本
        start_time = time.time()
        result_auto = random_walk_auto(P, rp=0.5, tol=0.01, use_mps=True)
        auto_time = time.time() - start_time
        
        acceleration_used = "MPS" if MPS_AVAILABLE else "CPU"
        print(f"✅ 自动版本完成（使用{acceleration_used}），耗时: {auto_time:.3f} 秒")
        
        if MPS_AVAILABLE and auto_time < cpu_time:
            speedup = cpu_time / auto_time
            print(f"🚀 MPS 加速比: {speedup:.2f}x")
        
        return True
        
    except Exception as e:
        print(f"❌ 插补函数测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """主测试函数"""
    print("MPS GPU 加速插补功能测试\n")
    
    # 测试 MPS 可用性
    mps_available = test_mps_availability()
    
    # 测试插补函数
    functions_work = test_impute_functions()
    
    print("\n=== 测试总结 ===")
    if mps_available:
        print("✅ MPS GPU 加速可用")
    else:
        print("❌ MPS GPU 加速不可用，将使用 CPU")
    
    if functions_work:
        print("✅ 插补函数正常工作")
    else:
        print("❌ 插补函数存在问题")
    
    if mps_available and functions_work:
        print("\n🎉 MPS GPU 加速插补功能已就绪！")
    else:
        print("\n⚠️  将使用 CPU 版本插补")

if __name__ == "__main__":
    main()
