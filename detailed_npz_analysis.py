#!/usr/bin/env python3
"""
详细分析 npz 文件内容并进行对比
"""

import numpy as np
import scipy.sparse as sp
from pathlib import Path
import matplotlib.pyplot as plt

def load_sparse_matrix_from_npz(filepath):
    """从 npz 文件加载稀疏矩阵"""
    data = np.load(filepath, allow_pickle=True)
    
    format_type = data['format'].item().decode('utf-8')
    shape = tuple(data['shape'])
    matrix_data = data['data']
    
    if format_type == 'coo':
        row = data['row']
        col = data['col']
        matrix = sp.coo_matrix((matrix_data, (row, col)), shape=shape)
    elif format_type == 'csr':
        indices = data['indices']
        indptr = data['indptr']
        matrix = sp.csr_matrix((matrix_data, indices, indptr), shape=shape)
    else:
        raise ValueError(f"Unsupported format: {format_type}")
    
    data.close()
    return matrix

def analyze_matrix_detailed(matrix, filename):
    """详细分析矩阵内容"""
    print(f"\n{'='*60}")
    print(f"详细分析: {filename}")
    print(f"{'='*60}")
    
    print(f"矩阵格式: {type(matrix).__name__}")
    print(f"矩阵形状: {matrix.shape}")
    print(f"数据类型: {matrix.dtype}")
    print(f"非零元素数: {matrix.nnz:,}")
    print(f"总元素数: {matrix.shape[0] * matrix.shape[1]:,}")
    print(f"稀疏度: {(1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1])) * 100:.2f}%")
    
    # 转换为密集矩阵的小样本进行分析
    if matrix.shape[0] > 0 and matrix.shape[1] > 0:
        sample_size = min(10, matrix.shape[0], matrix.shape[1])
        # 对于 COO 格式，先转换为 CSR 格式再切片
        if sp.isspmatrix_coo(matrix):
            matrix_csr = matrix.tocsr()
            sample = matrix_csr[:sample_size, :sample_size].toarray()
        else:
            sample = matrix[:sample_size, :sample_size].toarray()
        print(f"\n左上角 {sample_size}x{sample_size} 样本:")
        print(sample)
        
        # 数值统计
        if matrix.nnz > 0:
            dense_data = matrix.data
            print(f"\n数值统计:")
            print(f"  最小值: {np.min(dense_data):.6f}")
            print(f"  最大值: {np.max(dense_data):.6f}")
            print(f"  均值: {np.mean(dense_data):.6f}")
            print(f"  中位数: {np.median(dense_data):.6f}")
            print(f"  标准差: {np.std(dense_data):.6f}")
            
            # 数值分布
            unique_vals, counts = np.unique(dense_data, return_counts=True)
            print(f"  唯一值数量: {len(unique_vals)}")
            if len(unique_vals) <= 20:
                print(f"  所有唯一值: {unique_vals}")
                print(f"  对应计数: {counts}")
            else:
                print(f"  前10个最常见值:")
                sorted_indices = np.argsort(counts)[::-1][:10]
                for i in sorted_indices:
                    print(f"    值 {unique_vals[i]:.6f}: {counts[i]} 次")
    
    # 检查矩阵是否对称
    if matrix.shape[0] == matrix.shape[1]:
        # 对于大矩阵，只检查一个小的子矩阵
        check_size = min(100, matrix.shape[0])
        if sp.isspmatrix_coo(matrix):
            sub_matrix = matrix.tocsr()[:check_size, :check_size]
        else:
            sub_matrix = matrix[:check_size, :check_size]
        is_symmetric = np.allclose(sub_matrix.toarray(), sub_matrix.T.toarray())
        print(f"\n矩阵对称性检查 (前{check_size}x{check_size}): {'对称' if is_symmetric else '不对称'}")
    
    return matrix

def main():
    data_dir = Path("/Volumes/SumSung500/CSU/0_HiRES/data/example_data")
    npz_files = list(data_dir.glob("*.npz"))
    
    matrices = {}
    
    print("加载并分析所有矩阵...")
    for npz_file in npz_files:
        try:
            matrix = load_sparse_matrix_from_npz(str(npz_file))
            matrices[npz_file.name] = matrix
            analyze_matrix_detailed(matrix, npz_file.name)
        except Exception as e:
            print(f"Error loading {npz_file.name}: {e}")
    
    # 文件对比总结
    print(f"\n\n{'='*80}")
    print("文件对比总结")
    print(f"{'='*80}")
    
    print("\n1. 基本信息对比:")
    print(f"{'文件名':<50} {'形状':<15} {'非零元素':<12} {'稀疏度%':<8} {'数据类型':<10}")
    print("-" * 95)
    
    for name, matrix in matrices.items():
        sparsity = (1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1])) * 100
        print(f"{name:<50} {str(matrix.shape):<15} {matrix.nnz:<12,} {sparsity:<8.2f} {str(matrix.dtype):<10}")
    
    print("\n2. 数据特征对比:")
    for name, matrix in matrices.items():
        if matrix.nnz > 0:
            print(f"\n{name}:")
            data = matrix.data
            print(f"  数值范围: [{np.min(data):.6f}, {np.max(data):.6f}]")
            print(f"  均值±标准差: {np.mean(data):.6f}±{np.std(data):.6f}")
            print(f"  零值比例: {np.sum(data == 0) / len(data) * 100:.2f}%")
    
    print("\n3. 文件用途推测:")
    for name, matrix in matrices.items():
        print(f"\n{name}:")
        if "Astro_Endo_ODC_OPC" in name:
            print("  -> 可能是多种细胞类型(星形胶质细胞、内皮细胞、ODC、OPC)的染色体接触矩阵")
            print("  -> upper_tri 表示上三角矩阵，通常用于存储对称矩阵的一半")
        elif "GasdE701001" in name:
            print("  -> 可能是单细胞样本的染色体接触矩阵")
            print("  -> pad1_std1.0_rp0.5_sqrtvc 表示经过填充、标准化、随机化和方差归一化处理")
        elif "K562" in name:
            print("  -> K562 是常用的白血病细胞系")
            if "sim" in name:
                print("  -> 'sim' 可能表示模拟/仿真数据")
            elif "true" in name:
                print("  -> 'true' 可能表示真实/观测数据")
            print("  -> 2k 可能表示2kb分辨率")

if __name__ == "__main__":
    main()
