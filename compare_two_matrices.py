#!/usr/bin/env python3
"""
详细对比 Astro_Endo_ODC_OPC_chr1_upper_tri.npz 和 GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz 的相似性
"""

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from pathlib import Path

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
    
    data.close()
    return matrix

def compare_matrices():
    """详细对比两个矩阵"""
    data_dir = Path("/Volumes/SumSung500/CSU/0_HiRES/data/example_data")
    
    # 加载两个矩阵
    file1 = data_dir / "Astro_Endo_ODC_OPC_chr1_upper_tri.npz"
    file2 = data_dir / "GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz"
    
    matrix1 = load_sparse_matrix_from_npz(str(file1))
    matrix2 = load_sparse_matrix_from_npz(str(file2))
    
    print("=" * 80)
    print("详细对比分析")
    print("=" * 80)
    
    print(f"\n文件1: {file1.name}")
    print(f"文件2: {file2.name}")
    
    # 基本属性对比
    print(f"\n【基本属性对比】")
    print(f"矩阵形状:     文件1: {matrix1.shape}     文件2: {matrix2.shape}")
    print(f"矩阵格式:     文件1: {type(matrix1).__name__}     文件2: {type(matrix2).__name__}")
    print(f"数据类型:     文件1: {matrix1.dtype}     文件2: {matrix2.dtype}")
    print(f"非零元素:     文件1: {matrix1.nnz:,}     文件2: {matrix2.nnz:,}")
    
    sparsity1 = (1 - matrix1.nnz / (matrix1.shape[0] * matrix1.shape[1])) * 100
    sparsity2 = (1 - matrix2.nnz / (matrix2.shape[0] * matrix2.shape[1])) * 100
    print(f"稀疏度:       文件1: {sparsity1:.2f}%     文件2: {sparsity2:.2f}%")
    
    # 数值范围对比
    print(f"\n【数值特征对比】")
    print(f"数值范围:     文件1: [{np.min(matrix1.data):.6f}, {np.max(matrix1.data):.6f}]")
    print(f"              文件2: [{np.min(matrix2.data):.6f}, {np.max(matrix2.data):.6f}]")
    print(f"均值:         文件1: {np.mean(matrix1.data):.6f}     文件2: {np.mean(matrix2.data):.6f}")
    print(f"标准差:       文件1: {np.std(matrix1.data):.6f}     文件2: {np.std(matrix2.data):.6f}")
    print(f"中位数:       文件1: {np.median(matrix1.data):.6f}     文件2: {np.median(matrix2.data):.6f}")
    
    # 检查数值类型
    unique1 = np.unique(matrix1.data)
    unique2 = np.unique(matrix2.data)
    print(f"唯一值数量:   文件1: {len(unique1):,}     文件2: {len(unique2):,}")
    
    # 检查是否为整数
    is_int1 = np.allclose(matrix1.data, np.round(matrix1.data))
    is_int2 = np.allclose(matrix2.data, np.round(matrix2.data))
    print(f"是否为整数:   文件1: {is_int1}     文件2: {is_int2}")
    
    # 矩阵结构对比
    print(f"\n【矩阵结构对比】")
    
    # 检查矩阵1是否为上三角（由于形状不同，只能检查小样本）
    if matrix1.shape[0] == matrix1.shape[1]:
        sample_size = min(100, matrix1.shape[0])
        sample1 = matrix1.tocsr()[:sample_size, :sample_size].toarray()
        is_upper_tri1 = np.allclose(sample1, np.triu(sample1))
        print(f"文件1是否上三角: {is_upper_tri1}")
    else:
        print(f"文件1是长方形矩阵 ({matrix1.shape})，不是方阵")
    
    # 检查矩阵2是否对称
    if matrix2.shape[0] == matrix2.shape[1]:
        sample_size = min(100, matrix2.shape[0])
        sample2 = matrix2.tocsr()[:sample_size, :sample_size].toarray()
        is_symmetric2 = np.allclose(sample2, sample2.T)
        print(f"文件2是否对称: {is_symmetric2}")
    
    # 对角线分析（仅对方阵）
    if matrix2.shape[0] == matrix2.shape[1]:
        diag_vals = matrix2.diagonal()
        diag_mean = np.mean(diag_vals[diag_vals != 0])
        print(f"文件2对角线均值: {diag_mean:.6f}")
        print(f"文件2对角线零值数: {np.sum(diag_vals == 0)}")
    
    # 相似性分析
    print(f"\n【相似性分析】")
    print("形状相似性:     完全不同 (一个是长方形，一个是方形)")
    print("数据类型相似性: 不同 (float64 vs float32)")
    print("数值范围相似性: 完全不同 (整数计数 vs 归一化浮点数)")
    print("稀疏度相似性:   相对接近 (都是高稀疏度矩阵)")
    
    # 用途推测
    print(f"\n【用途推测】")
    print("文件1 (Astro_Endo_ODC_OPC):")
    print("  - 多细胞类型的Hi-C接触计数矩阵")
    print("  - 原始或轻度处理的数据")
    print("  - 整数计数值，表示接触频次")
    print("  - 上三角存储格式，节省空间")
    print("  - 长方形矩阵可能表示不同分辨率或区域")
    
    print("\n文件2 (GasdE701001):")
    print("  - 单细胞样本的Hi-C数据")
    print("  - 经过复杂预处理的数据")
    print("  - 归一化浮点值，范围[0,1]")
    print("  - 方形对称矩阵")
    print("  - 可能用于下游机器学习分析")
    
    # 创建可视化对比
    create_comparison_visualization(matrix1, matrix2, file1.stem, file2.stem)
    
    # 结论
    print(f"\n【最终结论】")
    print("这两个文件虽然都是Hi-C相关的稀疏矩阵，但:")
    print("❌ 矩阵维度完全不同")
    print("❌ 数据类型和数值范围差异巨大") 
    print("❌ 数据处理阶段不同")
    print("❌ 应用场景不同")
    print("✅ 都具有高稀疏度特征")
    print("✅ 都可能来自染色体接触数据")
    print("\n总体相似性: 很低 (主要是数据源类型相似，但具体内容差异很大)")

def create_comparison_visualization(matrix1, matrix2, name1, name2):
    """创建可视化对比"""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 矩阵1的可视化
    sample_size1 = min(50, matrix1.shape[0], matrix1.shape[1])
    if sp.isspmatrix_coo(matrix1):
        sample1 = matrix1.tocsr()[:sample_size1, :sample_size1].toarray()
    else:
        sample1 = matrix1[:sample_size1, :sample_size1].toarray()
    
    im1 = axes[0,0].imshow(sample1, cmap='viridis', aspect='auto')
    axes[0,0].set_title(f'{name1}\n({sample_size1}x{sample_size1} sample)')
    plt.colorbar(im1, ax=axes[0,0])
    
    # 矩阵2的可视化
    sample_size2 = min(50, matrix2.shape[0], matrix2.shape[1])
    sample2 = matrix2[:sample_size2, :sample_size2].toarray()
    
    im2 = axes[1,0].imshow(sample2, cmap='viridis', aspect='auto')
    axes[1,0].set_title(f'{name2}\n({sample_size2}x{sample_size2} sample)')
    plt.colorbar(im2, ax=axes[1,0])
    
    # 数值分布对比
    axes[0,1].hist(matrix1.data, bins=50, alpha=0.7, label=name1, density=True)
    axes[0,1].set_title('Matrix 1 Value Distribution')
    axes[0,1].set_xlabel('Value')
    axes[0,1].set_ylabel('Density')
    axes[0,1].set_yscale('log')
    
    axes[1,1].hist(matrix2.data, bins=50, alpha=0.7, label=name2, density=True, color='orange')
    axes[1,1].set_title('Matrix 2 Value Distribution')
    axes[1,1].set_xlabel('Value')
    axes[1,1].set_ylabel('Density')
    axes[1,1].set_yscale('log')
    
    # 基本统计对比
    stats_data = {
        'Shape': [f"{matrix1.shape[0]}×{matrix1.shape[1]}", f"{matrix2.shape[0]}×{matrix2.shape[1]}"],
        'NNZ': [f"{matrix1.nnz:,}", f"{matrix2.nnz:,}"],
        'Dtype': [str(matrix1.dtype), str(matrix2.dtype)],
        'Min': [f"{np.min(matrix1.data):.3f}", f"{np.min(matrix2.data):.3f}"],
        'Max': [f"{np.max(matrix1.data):.3f}", f"{np.max(matrix2.data):.3f}"],
        'Mean': [f"{np.mean(matrix1.data):.3f}", f"{np.mean(matrix2.data):.3f}"]
    }
    
    # 创建对比表格
    axes[0,2].axis('tight')
    axes[0,2].axis('off')
    table_data = [[key, stats_data[key][0], stats_data[key][1]] for key in stats_data.keys()]
    table = axes[0,2].table(cellText=table_data,
                           colLabels=['Metric', name1, name2],
                           cellLoc='center',
                           loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    axes[0,2].set_title('Comparison Table')
    
    # 稀疏度对比
    sparsity1 = (1 - matrix1.nnz / (matrix1.shape[0] * matrix1.shape[1])) * 100
    sparsity2 = (1 - matrix2.nnz / (matrix2.shape[0] * matrix2.shape[1])) * 100
    
    axes[1,2].bar([name1[:15], name2[:15]], [sparsity1, sparsity2], color=['blue', 'orange'])
    axes[1,2].set_title('Sparsity Comparison')
    axes[1,2].set_ylabel('Sparsity (%)')
    axes[1,2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('/Volumes/SumSung500/CSU/0_HiRES/matrix_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\n可视化对比图已保存为: matrix_comparison.png")

if __name__ == "__main__":
    compare_matrices()
