#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h5ad文件生成逻辑实现 - 修复版本（无psutil依赖）
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import gzip
import re
from collections import defaultdict
from scipy.sparse import csr_matrix, vstack
import gc


def parse_pairs_file(filepath, bin_size=20000):
    """
    解析.pairs.gz文件，生成接触矩阵
    
    Args:
        filepath (str): .pairs.gz文件路径
        bin_size (int): 分箱大小，默认20kb
        
    Returns:
        dict: 包含接触数据的字典
    """
    contacts = defaultdict(int)
    
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
                
            # 解析染色体和位置信息
            chrom1, pos1, chrom2, pos2 = parts[1], int(parts[2]), parts[3], int(parts[4])
            
            # 只处理相同染色体的接触
            if chrom1 == chrom2:
                # 计算bin索引
                bin1 = pos1 // bin_size
                bin2 = pos2 // bin_size
                
                # 确保bin1 <= bin2 (上三角矩阵)
                if bin1 > bin2:
                    bin1, bin2 = bin2, bin1
                    
                # 累加接触计数
                contacts[(bin1, bin2)] += 1
    
    return contacts


def create_feature_vector(contacts, max_bins=2000):
    """
    从接触数据创建特征向量（上三角矩阵部分）
    
    Args:
        contacts (dict): 接触数据字典
        max_bins (int): 最大bin数
        
    Returns:
        csr_matrix: 1×N的稀疏特征向量
    """
    # 计算上三角矩阵元素数量
    n_elements = max_bins * (max_bins + 1) // 2
    rows, cols, data = [], [], []
    
    for (bin1, bin2), count in contacts.items():
        if bin1 < max_bins and bin2 < max_bins and bin1 <= bin2:
            # 上三角矩阵索引映射
            # idx = bin1 * max_bins - bin1 * (bin1 - 1) // 2 + (bin2 - bin1)
            idx = int(bin1 * max_bins - bin1 * (bin1 - 1) // 2 + (bin2 - bin1))
            if idx < n_elements:
                rows.append(0)  # 只有一行
                cols.append(idx)
                data.append(count)
    
    return csr_matrix((data, (rows, cols)), shape=(1, n_elements))


def extract_cellname_from_filename(filename):
    """
    从scHi-C文件名中提取Cellname
    
    Args:
        filename (str): scHi-C文件名
        
    Returns:
        str: 提取的Cellname
    """
    match = re.search(r'_([^_]+)\.pairs\.gz$', filename)
    if match:
        return match.group(1)
    return None


def process_stage_to_h5ad(stage_files, metadata_df, output_file, max_bins=2000):
    """
    将一个Stage的数据处理为h5ad文件
    
    Args:
        stage_files (list): 该Stage的文件信息列表
        metadata_df (DataFrame): metadata数据框
        output_file (str): 输出h5ad文件路径
        max_bins (int): 最大bin数
    """
    print(f"开始处理Stage，共{len(stage_files)}个细胞文件")
    print(f"使用bin数: {max_bins}")
    
    # 收集所有细胞的特征向量
    all_vectors = []
    cell_names = []
    cell_metadata = []
    
    # 处理每个细胞文件
    for i, file_info in enumerate(stage_files):
        filename = file_info['filename']
        filepath = file_info['filepath']
        cellname = file_info['cellname']
        
        if i % 50 == 0:
            print(f"  处理进度: {i}/{len(stage_files)} - {filename}")
        
        try:
            # 解析pairs文件
            contacts = parse_pairs_file(filepath)
            
            # 创建特征向量
            vector = create_feature_vector(contacts, max_bins)
            
            # 保存向量和元数据
            all_vectors.append(vector)
            cell_names.append(cellname)
            
            # 获取细胞的metadata
            metadata_row = metadata_df[metadata_df['Cellname'] == cellname].iloc[0]
            cell_metadata.append(metadata_row.to_dict())
            
            # 定期清理内存
            if i % 200 == 0:
                gc.collect()
            
        except Exception as e:
            print(f"  处理文件 {filename} 时出错: {e}")
            continue
    
    if not all_vectors:
        print("  没有成功处理任何细胞文件")
        return
    
    print(f"  成功处理 {len(all_vectors)} 个细胞")
    
    # 合并所有特征向量
    try:
        print("  开始合并特征向量...")
        combined_matrix = vstack(all_vectors)
        print(f"  合并完成，矩阵形状: {combined_matrix.shape}")
    except Exception as e:
        print(f"  合并特征向量时出错: {e}")
        return
    
    # 创建obs数据框
    obs_df = pd.DataFrame(cell_metadata)
    obs_df.index = cell_names
    
    # 创建var数据框 (bin信息)
    # 使用上三角矩阵的索引命名
    n_elements = max_bins * (max_bins + 1) // 2
    var_df = pd.DataFrame(index=[f"bin_pair_{i}" for i in range(min(n_elements, combined_matrix.shape[1]))])
    
    # 创建AnnData对象
    print("  创建AnnData对象...")
    try:
        adata = sc.AnnData(X=combined_matrix, obs=obs_df, var=var_df)
        print(f"  AnnData对象创建完成，形状: {adata.shape}")
    except Exception as e:
        print(f"  创建AnnData对象时出错: {e}")
        return
    
    # 保存到h5ad文件
    try:
        print(f"  保存到文件: {output_file}")
        adata.write(output_file)
        print(f"  Stage数据已保存到: {output_file}")
    except Exception as e:
        print(f"  保存文件时出错: {e}")
    
    # 清理内存
    del all_vectors, combined_matrix, adata
    gc.collect()
    print(f"  处理完成")


def process_all_stages(stage_files_mapping, metadata_file, output_dir):
    """
    处理所有Stage的数据
    
    Args:
        stage_files_mapping (str): Stage文件映射CSV文件路径
        metadata_file (str): metadata文件路径
        output_dir (str): 输出目录路径
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取stage文件映射
    stage_df = pd.read_csv(stage_files_mapping)
    
    # 读取metadata
    metadata_df = pd.read_excel(metadata_file)
    
    # 按Stage分组处理
    for stage in stage_df['stage'].unique():
        print(f"\n{'='*50}")
        print(f"处理Stage: {stage}")
        print(f"{'='*50}")
        
        # 获取该Stage的所有文件
        stage_files = stage_df[stage_df['stage'] == stage].to_dict('records')
        
        # 定义输出文件路径
        output_file = os.path.join(output_dir, f"stage_{stage}.h5ad")
        
        # 处理该Stage的数据
        process_stage_to_h5ad(stage_files, metadata_df, output_file)


def main():
    """主函数"""
    stage_files_mapping = "/home/duxuyan/Projects/0_HiRES/code/claude_new/stage_files_mapping.csv"
    metadata_file = "/home/duxuyan/Projects/0_HiRES/data/GSE223917_HiRES_emb_metadata.xlsx"
    output_dir = "/home/duxuyan/Projects/0_HiRES/output"
    
    # 处理所有Stage的数据
    process_all_stages(stage_files_mapping, metadata_file, output_dir)


if __name__ == "__main__":
    main()