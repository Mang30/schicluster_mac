#!/usr/bin/env python3
"""
检查H5AD文件内容的脚本
"""

import anndata as ad
import pandas as pd
import numpy as np
import os

def check_h5ad_file(file_path):
    """检查单个H5AD文件的内容"""
    print(f"=== 检查文件: {os.path.basename(file_path)} ===")
    
    try:
        # 加载h5ad文件
        adata = ad.read_h5ad(file_path)
        
        print(f"AnnData对象信息:")
        print(f"  观察值 (细胞) 数量: {adata.n_obs}")
        print(f"  变量 (特征) 数量: {adata.n_vars}")
        print(f"  数据矩阵形状: {adata.X.shape}")
        print(f"  数据矩阵类型: {type(adata.X)}")
        
        if hasattr(adata.X, 'dtype'):
            print(f"  数据类型: {adata.X.dtype}")
        
        # 检查稀疏矩阵信息
        if hasattr(adata.X, 'nnz'):
            total_elements = adata.X.shape[0] * adata.X.shape[1]
            sparsity = adata.X.nnz / total_elements * 100
            print(f"  非零元素数量: {adata.X.nnz:,}")
            print(f"  稀疏度: {sparsity:.6f}%")
        
        # 检查观察值信息 (细胞元数据)
        print(f"\n--- 细胞元数据 (obs) ---")
        print(f"  obs 列数: {adata.obs.shape[1]}")
        if adata.obs.shape[1] > 0:
            print(f"  obs 列名: {list(adata.obs.columns)}")
            print(f"  示例数据:")
            print(adata.obs.head(3))
        
        # 检查变量信息 (特征元数据)
        print(f"\n--- 特征元数据 (var) ---")
        print(f"  var 列数: {adata.var.shape[1]}")
        if adata.var.shape[1] > 0:
            print(f"  var 列名: {list(adata.var.columns)}")
            print(f"  示例数据:")
            print(adata.var.head(3))
        
        # 检查数据矩阵的一小部分
        print(f"\n--- 数据矩阵示例 ---")
        if hasattr(adata.X, 'toarray'):
            # 稀疏矩阵
            sample_data = adata.X[:min(5, adata.n_obs), :min(10, adata.n_vars)].toarray()
        else:
            # 密集矩阵
            sample_data = adata.X[:min(5, adata.n_obs), :min(10, adata.n_vars)]
        
        print("左上角5x10区域:")
        print(sample_data)
        
        # 检查数据统计
        if hasattr(adata.X, 'data'):
            # 稀疏矩阵
            data_values = adata.X.data
            print(f"\n--- 数据统计 ---")
            print(f"非零值范围: [{data_values.min():.6f}, {data_values.max():.6f}]")
            print(f"非零值均值: {data_values.mean():.6f}")
            print(f"非零值标准差: {data_values.std():.6f}")
        
        # 检查是否有其他存储的层
        if adata.layers:
            print(f"\n--- 数据层 (layers) ---")
            for layer_name in adata.layers.keys():
                layer = adata.layers[layer_name]
                print(f"  {layer_name}: {layer.shape}, {type(layer)}")
        
        # 检查是否有无结构数据
        if adata.uns:
            print(f"\n--- 无结构数据 (uns) ---")
            for key in adata.uns.keys():
                value = adata.uns[key]
                print(f"  {key}: {type(value)}")
                if isinstance(value, (str, int, float)):
                    print(f"    值: {value}")
                elif hasattr(value, 'shape'):
                    print(f"    形状: {value.shape}")
        
        return True
        
    except Exception as e:
        print(f"❌ 加载文件出错: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """主函数"""
    # 检查Hires_h5ad_data目录下的文件
    data_dir = "/Volumes/SumSung500/CSU/0_HiRES/data/example_data/Hires_h5ad_data"
    
    print("🔍 搜索H5AD文件...")
    
    if not os.path.exists(data_dir):
        print(f"❌ 数据目录不存在: {data_dir}")
        return
    
    # 获取所有h5ad文件
    h5ad_files = [f for f in os.listdir(data_dir) if f.endswith('.h5ad')]
    
    if not h5ad_files:
        print("❌ 没有找到任何H5AD文件")
        return
    
    print(f"✅ 找到 {len(h5ad_files)} 个H5AD文件")
    
    # 检查前几个文件作为示例
    max_files_to_check = 2
    for i, filename in enumerate(sorted(h5ad_files)[:max_files_to_check]):
        file_path = os.path.join(data_dir, filename)
        print(f"\n{'='*60}")
        print(f"文件 {i+1}/{min(len(h5ad_files), max_files_to_check)}")
        success = check_h5ad_file(file_path)
        
        if not success:
            continue
    
    if len(h5ad_files) > max_files_to_check:
        print(f"\n... 还有 {len(h5ad_files) - max_files_to_check} 个文件未显示")

if __name__ == "__main__":
    main()
