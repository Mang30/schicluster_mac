#!/usr/bin/env python3
"""
分析 data/example_data/ 目录下的 npz 文件内容
"""

import numpy as np
import os
import scipy.sparse as sp
from pathlib import Path

def analyze_npz_file(filepath):
    """分析单个 npz 文件"""
    print(f"\n{'='*80}")
    print(f"文件: {os.path.basename(filepath)}")
    print(f"路径: {filepath}")
    print(f"{'='*80}")
    
    try:
        # 加载 npz 文件
        data = np.load(filepath, allow_pickle=True)
        
        print(f"文件大小: {os.path.getsize(filepath) / (1024*1024):.2f} MB")
        print(f"包含的数组/对象数量: {len(data.files)}")
        print(f"数组名称: {list(data.files)}")
        
        # 分析每个数组
        for key in data.files:
            print(f"\n--- 数组: '{key}' ---")
            arr = data[key]
            
            # 检查是否是稀疏矩阵
            if hasattr(arr, 'item') and callable(getattr(arr, 'item')):
                # 可能是包装的稀疏矩阵
                try:
                    arr_item = arr.item()
                    if sp.issparse(arr_item):
                        print(f"类型: 稀疏矩阵 ({type(arr_item).__name__})")
                        print(f"形状: {arr_item.shape}")
                        print(f"数据类型: {arr_item.dtype}")
                        print(f"非零元素数: {arr_item.nnz}")
                        print(f"稀疏度: {(1 - arr_item.nnz / (arr_item.shape[0] * arr_item.shape[1])) * 100:.2f}%")
                        if arr_item.shape[0] > 0 and arr_item.shape[1] > 0:
                            dense_sample = arr_item[:min(5, arr_item.shape[0]), :min(5, arr_item.shape[1])].toarray()
                            print(f"左上角 5x5 样本:\n{dense_sample}")
                    else:
                        print(f"类型: {type(arr_item)}")
                        print(f"内容: {arr_item}")
                except:
                    print(f"类型: {type(arr)}")
                    print(f"形状: {getattr(arr, 'shape', 'N/A')}")
                    print(f"数据类型: {getattr(arr, 'dtype', type(arr))}")
            else:
                print(f"类型: {type(arr)} (numpy数组)")
                print(f"形状: {arr.shape}")
                print(f"数据类型: {arr.dtype}")
                
                if arr.ndim <= 2 and arr.size <= 25:
                    print(f"内容:\n{arr}")
                elif arr.ndim == 1:
                    print(f"前10个元素: {arr[:10]}")
                    if len(arr) > 10:
                        print(f"后10个元素: {arr[-10:]}")
                elif arr.ndim == 2:
                    print(f"左上角样本 (最多5x5):\n{arr[:min(5, arr.shape[0]), :min(5, arr.shape[1])]}")
                else:
                    print(f"多维数组，形状: {arr.shape}")
                
                # 统计信息
                if np.issubdtype(arr.dtype, np.number):
                    print(f"数值统计:")
                    print(f"  最小值: {np.min(arr)}")
                    print(f"  最大值: {np.max(arr)}")
                    print(f"  均值: {np.mean(arr):.6f}")
                    print(f"  标准差: {np.std(arr):.6f}")
                    print(f"  零值数量: {np.sum(arr == 0)}")
                    print(f"  非零值数量: {np.sum(arr != 0)}")
        
        data.close()
        
    except Exception as e:
        print(f"Error loading {filepath}: {e}")

def main():
    # 分析目录下的所有 npz 文件
    data_dir = Path("/Volumes/SumSung500/CSU/0_HiRES/data/example_data")
    npz_files = list(data_dir.glob("*.npz"))
    
    print("发现以下 npz 文件:")
    for f in npz_files:
        print(f"  - {f.name}")
    
    # 逐个分析
    for npz_file in npz_files:
        analyze_npz_file(str(npz_file))
    
    # 总结对比
    print(f"\n\n{'='*80}")
    print("文件对比总结")
    print(f"{'='*80}")
    
    print(f"共发现 {len(npz_files)} 个 npz 文件:")
    for f in npz_files:
        size_mb = os.path.getsize(f) / (1024*1024)
        print(f"  - {f.name}: {size_mb:.2f} MB")

if __name__ == "__main__":
    main()
