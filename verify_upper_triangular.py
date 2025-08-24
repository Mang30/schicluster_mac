#!/usr/bin/env python3
"""
验证上三角矩阵转换结果
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

def verify_upper_triangular_conversion():
    """验证上三角矩阵转换结果"""
    data_dir = Path("/Volumes/SumSung500/CSU/0_HiRES/data/example_data")
    
    # 文件路径
    original_file = data_dir / "GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz"
    upper_tri_file = data_dir / "GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc_upper_tri.npz"
    upper_tri_no_diag_file = data_dir / "GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc_upper_tri_no_diag.npz"
    
    print("=" * 80)
    print("验证上三角矩阵转换结果")
    print("=" * 80)
    
    # 加载所有矩阵
    print("加载矩阵...")
    original = load_sparse_matrix_from_npz(str(original_file))
    upper_tri = load_sparse_matrix_from_npz(str(upper_tri_file))
    upper_tri_no_diag = load_sparse_matrix_from_npz(str(upper_tri_no_diag_file))
    
    print(f"\n文件信息:")
    print(f"原始矩阵:           {original.shape}, {original.nnz:,} 非零元素")
    print(f"上三角矩阵(含对角): {upper_tri.shape}, {upper_tri.nnz:,} 非零元素")
    print(f"上三角矩阵(无对角): {upper_tri_no_diag.shape}, {upper_tri_no_diag.nnz:,} 非零元素")
    
    # 验证上三角性质
    print(f"\n验证上三角性质...")
    
    # 转换为COO格式以便检查
    if not sp.isspmatrix_coo(upper_tri):
        upper_tri_coo = upper_tri.tocoo()
    else:
        upper_tri_coo = upper_tri
    
    if not sp.isspmatrix_coo(upper_tri_no_diag):
        upper_tri_no_diag_coo = upper_tri_no_diag.tocoo()
    else:
        upper_tri_no_diag_coo = upper_tri_no_diag
    
    # 检查上三角条件
    is_upper_tri = np.all(upper_tri_coo.row <= upper_tri_coo.col)
    is_upper_tri_no_diag = np.all(upper_tri_no_diag_coo.row < upper_tri_no_diag_coo.col)
    
    print(f"包含对角线版本是否满足上三角条件 (row <= col): {is_upper_tri}")
    print(f"不含对角线版本是否满足上三角条件 (row < col): {is_upper_tri_no_diag}")
    
    # 检查对角线元素
    diag_count_original = np.sum(upper_tri_coo.row == upper_tri_coo.col)
    diag_count_no_diag = np.sum(upper_tri_no_diag_coo.row == upper_tri_no_diag_coo.col)
    
    print(f"\n对角线元素数量:")
    print(f"包含对角线版本: {diag_count_original}")
    print(f"不含对角线版本: {diag_count_no_diag}")
    
    # 验证数据完整性 - 检查是否正确提取了上三角部分
    print(f"\n验证数据完整性...")
    
    # 从原始矩阵重构上三角部分进行对比
    original_coo = original.tocoo()
    
    # 原始矩阵的上三角部分（包含对角线）
    mask_with_diag = original_coo.row <= original_coo.col
    expected_nnz_with_diag = np.sum(mask_with_diag)
    
    # 原始矩阵的上三角部分（不含对角线）
    mask_no_diag = original_coo.row < original_coo.col
    expected_nnz_no_diag = np.sum(mask_no_diag)
    
    print(f"预期非零元素数 (含对角): {expected_nnz_with_diag:,}")
    print(f"实际非零元素数 (含对角): {upper_tri.nnz:,}")
    print(f"匹配: {expected_nnz_with_diag == upper_tri.nnz}")
    
    print(f"预期非零元素数 (无对角): {expected_nnz_no_diag:,}")
    print(f"实际非零元素数 (无对角): {upper_tri_no_diag.nnz:,}")
    print(f"匹配: {expected_nnz_no_diag == upper_tri_no_diag.nnz}")
    
    # 可视化对比
    create_visualization_comparison(original, upper_tri, upper_tri_no_diag)
    
    # 数据重构测试
    print(f"\n数据重构测试...")
    reconstructed = upper_tri + upper_tri.T
    
    # 需要减去对角线重复计算的部分
    diag_matrix = sp.diags(upper_tri.diagonal(), shape=upper_tri.shape, format='coo')
    reconstructed = reconstructed - diag_matrix
    
    # 检查重构精度
    sample_size = min(100, original.shape[0])
    orig_sample = original.tocsr()[:sample_size, :sample_size].toarray()
    recon_sample = reconstructed.tocsr()[:sample_size, :sample_size].toarray()
    
    reconstruction_error = np.mean(np.abs(orig_sample - recon_sample))
    print(f"重构误差 (前{sample_size}x{sample_size}): {reconstruction_error:.10f}")
    print(f"重构成功: {reconstruction_error < 1e-10}")

def create_visualization_comparison(original, upper_tri, upper_tri_no_diag):
    """创建可视化对比"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    sample_size = min(100, original.shape[0])
    
    # 原始矩阵
    orig_sample = original.tocsr()[:sample_size, :sample_size].toarray()
    im1 = axes[0].imshow(orig_sample, cmap='viridis', aspect='auto')
    axes[0].set_title(f'Original Matrix\n({sample_size}x{sample_size} sample)')
    plt.colorbar(im1, ax=axes[0])
    
    # 上三角矩阵（含对角线）
    upper_sample = upper_tri.tocsr()[:sample_size, :sample_size].toarray()
    im2 = axes[1].imshow(upper_sample, cmap='viridis', aspect='auto')
    axes[1].set_title(f'Upper Triangular (with diagonal)\n({sample_size}x{sample_size} sample)')
    plt.colorbar(im2, ax=axes[1])
    
    # 上三角矩阵（无对角线）
    upper_no_diag_sample = upper_tri_no_diag.tocsr()[:sample_size, :sample_size].toarray()
    im3 = axes[2].imshow(upper_no_diag_sample, cmap='viridis', aspect='auto')
    axes[2].set_title(f'Upper Triangular (no diagonal)\n({sample_size}x{sample_size} sample)')
    plt.colorbar(im3, ax=axes[2])
    
    plt.tight_layout()
    plt.savefig('/Volumes/SumSung500/CSU/0_HiRES/upper_triangular_verification.png', 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\n可视化对比图已保存为: upper_triangular_verification.png")

if __name__ == "__main__":
    verify_upper_triangular_conversion()
