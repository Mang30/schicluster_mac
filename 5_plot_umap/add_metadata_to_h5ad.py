#!/usr/bin/env python3
"""
1_add_metadata_to_h5ad.py - 为h5ad文件添加元数据
将5_plot_umap/GSE223917_HiRES_emb_metadata.xlsx中的元数据添加到h5ad文件中
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import re

def extract_cellname_from_obsname(obs_name):
    """
    从观测名称中提取Cellname
    例如: GSM7000434_GaseE752203 -> GaseE752203
    """
    # 使用正则表达式提取最后的Cellname部分
    match = re.search(r'([^_]+)$', obs_name)
    if match:
        return match.group(1)
    return obs_name

def add_metadata_to_h5ad(h5ad_file, metadata_file, output_dir):
    """
    为h5ad文件添加元数据
    """
    print(f"处理文件: {h5ad_file}")
    
    # 读取h5ad文件
    adata = sc.read(h5ad_file)
    print(f"原始h5ad形状: {adata.shape}")
    print(f"原始观测列: {list(adata.obs.columns)}")
    
    # 读取元数据
    metadata_df = pd.read_excel(metadata_file)
    print(f"元数据形状: {metadata_df.shape}")
    
    # 创建Cellname到元数据的映射
    metadata_dict = {}
    for _, row in metadata_df.iterrows():
        metadata_dict[row['Cellname']] = row
    
    # 为每个观测添加元数据
    new_obs_data = {}
    
    # 添加现有列
    for col in adata.obs.columns:
        new_obs_data[col] = adata.obs[col].tolist()
    
    # 添加新的元数据列
    metadata_columns = ['Stage', 'Celltype', 'Cellcycle phase', 'CDPS cluster', 'Sub_k_cluster']
    for col in metadata_columns:
        new_obs_data[col] = []
    
    # 匹配并添加元数据
    for obs_name in adata.obs_names:
        cellname = extract_cellname_from_obsname(obs_name)
        
        if cellname in metadata_dict:
            row = metadata_dict[cellname]
            for col in metadata_columns:
                new_obs_data[col].append(row[col])
        else:
            print(f"警告: 未找到细胞 {obs_name} (提取的Cellname: {cellname}) 的元数据")
            for col in metadata_columns:
                new_obs_data[col].append(np.nan)
    
    # 更新obs
    adata.obs = pd.DataFrame(new_obs_data, index=adata.obs_names)
    
    print(f"更新后的观测列: {list(adata.obs.columns)}")
    
    # 保存更新后的h5ad文件
    output_file = os.path.join(output_dir, os.path.basename(h5ad_file))
    adata.write(output_file)
    print(f"保存更新后的文件: {output_file}")
    
    return adata

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='为h5ad文件添加元数据')
    parser.add_argument('--input_dir', type=str, required=True, 
                        help='输入h5ad文件的目录路径')
    parser.add_argument('--metadata_file', type=str, required=True, 
                        help='元数据文件路径')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='输出h5ad文件的目录路径')
    parser.add_argument('--stage', type=str, required=False,
                        help='特定stage名称 (例如: stage_E75)。如果未指定，则处理所有stage文件')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 设置路径
    input_dir = args.input_dir
    metadata_file = args.metadata_file
    output_dir = args.output_dir
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 根据是否指定stage来确定要处理的文件
    if args.stage:
        # 处理特定stage的文件
        h5ad_files = list(Path(input_dir).glob(f"{args.stage}_decay_profiles.h5ad"))
        print(f"找到 {len(h5ad_files)} 个 {args.stage} 的h5ad文件")
    else:
        # 处理所有h5ad文件
        h5ad_files = list(Path(input_dir).glob("*.h5ad"))
        print(f"找到 {len(h5ad_files)} 个h5ad文件")
    
    for h5ad_file in h5ad_files:
        try:
            add_metadata_to_h5ad(str(h5ad_file), metadata_file, output_dir)
        except Exception as e:
            print(f"处理文件 {h5ad_file} 时出错: {e}")
    
    print("所有文件处理完成!")

if __name__ == "__main__":
    main()