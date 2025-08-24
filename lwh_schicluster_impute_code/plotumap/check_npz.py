#!/usr/bin/env python3
"""
检查NPZ文件内容的脚本
"""

import numpy as np
import os
from scipy.sparse import load_npz

def check_npz_file(file_path):
    """检查单个NPZ文件的内容"""
    print(f"=== 检查文件: {os.path.basename(file_path)} ===")
    
    try:
        # 方法1：直接用numpy加载查看内容结构
        with np.load(file_path) as data:
            print(f"NPZ文件包含的数组键: {list(data.keys())}")
            for key in data.keys():
                arr = data[key]
                print(f"  {key}: shape={arr.shape}, dtype={arr.dtype}")
                if hasattr(arr, 'size'):
                    print(f"    总元素数: {arr.size}")
        
        # 方法2：用scipy加载稀疏矩阵
        print(f"\n--- 稀疏矩阵详细信息 ---")
        matrix = load_npz(file_path)
        print(f"矩阵类型: {type(matrix)}")
        print(f"矩阵形状: {matrix.shape}")
        print(f"矩阵数据类型: {matrix.dtype}")
        print(f"非零元素数量: {matrix.nnz:,}")
        
        total_elements = matrix.shape[0] * matrix.shape[1]
        sparsity = matrix.nnz / total_elements * 100
        print(f"总元素数量: {total_elements:,}")
        print(f"稀疏度: {sparsity:.6f}%")
        print(f"零元素比例: {100-sparsity:.6f}%")
        
        # 查看数据范围
        if matrix.nnz > 0:
            print(f"\n--- 数据统计 ---")
            data = matrix.data  # 非零元素的值
            print(f"非零值范围: [{data.min():.3f}, {data.max():.3f}]")
            print(f"非零值均值: {data.mean():.3f}")
            print(f"非零值标准差: {data.std():.3f}")
        
        # 查看矩阵的一小部分数据示例
        print(f"\n--- 矩阵数据示例 (左上角10x10) ---")
        sample_size = min(10, matrix.shape[0], matrix.shape[1])
        dense_sample = matrix[:sample_size, :sample_size].toarray()
        print(dense_sample)
        
        return True
        
    except Exception as e:
        print(f"❌ 加载文件出错: {e}")
        return False

def main():
    """主函数"""
    # 查找已存在的merged NPZ文件
    output_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/merged_matrices"
    
    print("🔍 搜索NPZ文件...")
    
    if not os.path.exists(output_dir):
        print(f"❌ 输出目录不存在: {output_dir}")
        return
    
    # 递归搜索所有NPZ文件
    npz_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith('_merged.npz'):
                npz_files.append(os.path.join(root, file))
    
    if not npz_files:
        print("❌ 没有找到任何merged NPZ文件")
        return
    
    print(f"✅ 找到 {len(npz_files)} 个NPZ文件")
    
    # 检查前几个文件作为示例
    max_files_to_check = 3
    for i, file_path in enumerate(npz_files[:max_files_to_check]):
        print(f"\n{'='*60}")
        print(f"文件 {i+1}/{min(len(npz_files), max_files_to_check)}")
        success = check_npz_file(file_path)
        
        if not success:
            continue
    
    if len(npz_files) > max_files_to_check:
        print(f"\n... 还有 {len(npz_files) - max_files_to_check} 个文件未显示")

if __name__ == "__main__":
    main()
