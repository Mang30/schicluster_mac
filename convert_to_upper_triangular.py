#!/usr/bin/env python3
"""
将 GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz 中的方形矩阵转换为上三角矩阵
"""

import numpy as np
import scipy.sparse as sp
from pathlib import Path
import os

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

def save_sparse_matrix_to_npz(matrix, filepath, format_type='coo'):
    """将稀疏矩阵保存为 npz 文件"""
    if format_type == 'coo':
        if not sp.isspmatrix_coo(matrix):
            matrix = matrix.tocoo()
        
        np.savez_compressed(filepath,
                          row=matrix.row,
                          col=matrix.col,
                          data=matrix.data,
                          shape=np.array(matrix.shape),
                          format=format_type.encode('utf-8'))
    elif format_type == 'csr':
        if not sp.isspmatrix_csr(matrix):
            matrix = matrix.tocsr()
        
        np.savez_compressed(filepath,
                          indices=matrix.indices,
                          indptr=matrix.indptr,
                          data=matrix.data,
                          shape=np.array(matrix.shape),
                          format=format_type.encode('utf-8'))
    else:
        raise ValueError(f"Unsupported format: {format_type}")

def convert_to_upper_triangular(input_file, output_file=None, include_diagonal=True):
    """
    将方形矩阵转换为上三角矩阵
    
    Parameters:
    -----------
    input_file : str
        输入文件路径
    output_file : str, optional
        输出文件路径，如果为None则自动生成
    include_diagonal : bool, default=True
        是否包含对角线元素
    """
    
    print(f"正在加载矩阵: {input_file}")
    
    # 加载原始矩阵
    matrix = load_sparse_matrix_from_npz(input_file)
    
    print(f"原始矩阵信息:")
    print(f"  形状: {matrix.shape}")
    print(f"  格式: {type(matrix).__name__}")
    print(f"  数据类型: {matrix.dtype}")
    print(f"  非零元素: {matrix.nnz:,}")
    print(f"  稀疏度: {(1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1])) * 100:.2f}%")
    
    # 检查是否为方形矩阵
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"输入矩阵不是方形矩阵，形状为 {matrix.shape}")
    
    # 检查是否为对称矩阵（抽样检查）
    print(f"\n检查矩阵对称性...")
    sample_size = min(100, matrix.shape[0])
    if sp.isspmatrix_coo(matrix):
        sample_matrix = matrix.tocsr()[:sample_size, :sample_size]
    else:
        sample_matrix = matrix[:sample_size, :sample_size]
    
    is_symmetric = np.allclose(sample_matrix.toarray(), sample_matrix.T.toarray())
    print(f"矩阵是否对称 (前{sample_size}x{sample_size}): {is_symmetric}")
    
    # 转换为COO格式以便操作
    if not sp.isspmatrix_coo(matrix):
        print("转换为COO格式...")
        matrix = matrix.tocoo()
    
    print(f"\n开始提取上三角部分...")
    
    # 提取上三角部分
    if include_diagonal:
        # 包含对角线: row <= col
        mask = matrix.row <= matrix.col
        print("提取上三角矩阵（包含对角线）")
    else:
        # 不包含对角线: row < col
        mask = matrix.row < matrix.col
        print("提取上三角矩阵（不包含对角线）")
    
    # 应用掩码
    upper_tri_row = matrix.row[mask]
    upper_tri_col = matrix.col[mask]
    upper_tri_data = matrix.data[mask]
    
    # 创建新的上三角矩阵
    upper_tri_matrix = sp.coo_matrix((upper_tri_data, (upper_tri_row, upper_tri_col)), 
                                   shape=matrix.shape)
    
    print(f"\n上三角矩阵信息:")
    print(f"  形状: {upper_tri_matrix.shape}")
    print(f"  非零元素: {upper_tri_matrix.nnz:,}")
    print(f"  原矩阵非零元素: {matrix.nnz:,}")
    print(f"  保留比例: {upper_tri_matrix.nnz / matrix.nnz * 100:.2f}%")
    print(f"  新稀疏度: {(1 - upper_tri_matrix.nnz / (upper_tri_matrix.shape[0] * upper_tri_matrix.shape[1])) * 100:.2f}%")
    
    # 生成输出文件名
    if output_file is None:
        input_path = Path(input_file)
        base_name = input_path.stem
        if include_diagonal:
            output_file = input_path.parent / f"{base_name}_upper_tri.npz"
        else:
            output_file = input_path.parent / f"{base_name}_upper_tri_no_diag.npz"
    
    # 保存结果
    print(f"\n保存上三角矩阵到: {output_file}")
    save_sparse_matrix_to_npz(upper_tri_matrix, str(output_file), format_type='coo')
    
    # 验证保存的文件
    print(f"\n验证保存的文件...")
    verification_matrix = load_sparse_matrix_from_npz(str(output_file))
    print(f"验证成功: 形状 {verification_matrix.shape}, 非零元素 {verification_matrix.nnz:,}")
    
    # 计算文件大小变化
    original_size = os.path.getsize(input_file) / (1024 * 1024)  # MB
    new_size = os.path.getsize(output_file) / (1024 * 1024)  # MB
    print(f"\n文件大小变化:")
    print(f"  原文件: {original_size:.2f} MB")
    print(f"  新文件: {new_size:.2f} MB")
    print(f"  大小比例: {new_size / original_size * 100:.2f}%")
    
    return str(output_file)

def main():
    """主函数"""
    input_file = "/Volumes/SumSung500/CSU/0_HiRES/data/example_data/GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz"
    
    print("=" * 80)
    print("矩阵上三角转换工具")
    print("=" * 80)
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误: 输入文件不存在: {input_file}")
        return
    
    try:
        # 转换为上三角矩阵（包含对角线）
        output_file = convert_to_upper_triangular(input_file, include_diagonal=True)
        print(f"\n✅ 转换完成!")
        print(f"输出文件: {output_file}")
        
        # 也可以选择创建不包含对角线的版本
        print(f"\n" + "="*50)
        print("创建不包含对角线的版本...")
        output_file_no_diag = convert_to_upper_triangular(input_file, include_diagonal=False)
        print(f"\n✅ 无对角线版本创建完成!")
        print(f"输出文件: {output_file_no_diag}")
        
    except Exception as e:
        print(f"❌ 转换失败: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
