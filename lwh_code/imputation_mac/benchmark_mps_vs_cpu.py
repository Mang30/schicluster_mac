#!/usr/bin/env python3
"""
完整的 MPS vs CPU 性能基准测试
"""

import sys
import time
import logging
import numpy as np
from scipy.sparse import csr_matrix

# 设置日志级别
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def create_test_matrix(size, sparsity=0.1):
    """创建指定大小和稀疏度的测试矩阵"""
    np.random.seed(42)
    nnz = int(size * size * sparsity)
    data = np.random.random(nnz)
    row = np.random.randint(0, size, nnz)
    col = np.random.randint(0, size, nnz)
    P = csr_matrix((data, (row, col)), shape=(size, size))
    
    # 归一化为概率矩阵
    P = P + P.T  # 对称化
    row_sums = np.array(P.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1  # 避免除零
    P = P.multiply(1/row_sums[:, np.newaxis])
    
    return P

def benchmark_performance():
    """性能基准测试"""
    print("=== MPS vs CPU 性能基准测试 ===\n")
    
    try:
        # 导入函数
        sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')
        from schicluster.impute.impute_chromosome import (
            random_walk_cpu, 
            random_walk_mps,
            random_walk_auto,
            MPS_AVAILABLE
        )
        
        # 测试不同大小的矩阵
        test_sizes = [100, 500, 1000, 2000]
        results = []
        
        for size in test_sizes:
            print(f"📊 测试矩阵大小: {size}x{size}")
            
            # 创建测试矩阵
            P = create_test_matrix(size, sparsity=0.1)
            actual_sparsity = P.nnz / (size * size)
            print(f"   实际稀疏度: {actual_sparsity:.4f}")
            
            # CPU 测试
            start_time = time.time()
            result_cpu = random_walk_cpu(P, rp=0.5, tol=0.01)
            cpu_time = time.time() - start_time
            print(f"   CPU 时间: {cpu_time:.3f} 秒")
            
            # MPS 测试（如果可用）
            mps_time = None
            speedup = None
            if MPS_AVAILABLE:
                try:
                    start_time = time.time()
                    result_mps = random_walk_mps(P, rp=0.5, tol=0.01)
                    mps_time = time.time() - start_time
                    print(f"   MPS 时间: {mps_time:.3f} 秒")
                    
                    speedup = cpu_time / mps_time
                    if speedup > 1:
                        print(f"   🚀 MPS 加速: {speedup:.2f}x")
                    else:
                        print(f"   🐌 MPS 较慢: {1/speedup:.2f}x")
                except Exception as e:
                    print(f"   ❌ MPS 测试失败: {e}")
            
            # 自动选择测试
            start_time = time.time()
            result_auto = random_walk_auto(P, rp=0.5, tol=0.01, use_mps=True)
            auto_time = time.time() - start_time
            print(f"   自动选择时间: {auto_time:.3f} 秒")
            
            results.append({
                'size': size,
                'sparsity': actual_sparsity,
                'cpu_time': cpu_time,
                'mps_time': mps_time,
                'auto_time': auto_time,
                'speedup': speedup
            })
            print()
        
        # 总结报告
        print("=== 性能总结 ===")
        print(f"{'大小':<8} {'稀疏度':<8} {'CPU':<8} {'MPS':<8} {'自动':<8} {'加速比':<8}")
        print("-" * 60)
        
        for r in results:
            size_str = f"{r['size']}x{r['size']}"
            sparsity_str = f"{r['sparsity']:.3f}"
            cpu_str = f"{r['cpu_time']:.3f}s"
            mps_str = f"{r['mps_time']:.3f}s" if r['mps_time'] else "N/A"
            auto_str = f"{r['auto_time']:.3f}s"
            speedup_str = f"{r['speedup']:.2f}x" if r['speedup'] else "N/A"
            
            print(f"{size_str:<8} {sparsity_str:<8} {cpu_str:<8} {mps_str:<8} {auto_str:<8} {speedup_str:<8}")
        
        # 建议
        print("\n=== 使用建议 ===")
        best_mps_size = None
        for r in results:
            if r['speedup'] and r['speedup'] > 1:
                best_mps_size = r['size']
                break
        
        if best_mps_size:
            print(f"✅ 建议矩阵大小 >= {best_mps_size} 时使用 MPS 加速")
        else:
            print("⚠️  在测试的矩阵大小范围内，CPU 性能更佳")
            print("💡 尝试更大的矩阵（如 5000x5000）可能会看到 MPS 优势")
        
        return True
        
    except Exception as e:
        print(f"❌ 基准测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """主函数"""
    print("MPS vs CPU 完整性能基准测试\n")
    
    # 检查 MPS 可用性
    try:
        import torch
        print(f"PyTorch 版本: {torch.__version__}")
        if torch.backends.mps.is_available():
            print("✅ MPS 加速可用\n")
        else:
            print("❌ MPS 加速不可用\n")
    except ImportError:
        print("❌ PyTorch 未安装\n")
    
    # 运行基准测试
    benchmark_performance()

if __name__ == "__main__":
    main()
