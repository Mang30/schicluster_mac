#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h5ad文件生成逻辑实现
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import gzip
import re
from collections import defaultdict, Counter
from scipy.sparse import csr_matrix
import h5py


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


def create_sparse_matrix(contacts, max_bins=10000):
    """
    从接触数据创建稀疏矩阵
    
    Args:
        contacts (dict): 接触数据字典
        max_bins (int): 最大bin数
        
    Returns:
        csr_matrix: 稀疏接触矩阵
    """
    # 提取行、列和数据
    rows, cols, data = [], [], []
    
    for (bin1, bin2), count in contacts.items():
        if bin1 < max_bins and bin2 < max_bins:
            rows.append(bin1)
            cols.append(bin2)
            data.append(count)
            
            # 如果不是对角线元素，也要添加对称元素
            if bin1 != bin2:
                rows.append(bin2)
                cols.append(bin1)
                data.append(count)
    
    # 创建稀疏矩阵
    matrix = csr_matrix((data, (rows, cols)), shape=(max_bins, max_bins))
    return matrix


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


def process_stage_to_h5ad(stage_files, metadata_df, output_file, max_bins=10000):
    """
    将一个Stage的数据处理为h5ad文件
    
    Args:
        stage_files (list): 该Stage的文件信息列表
        metadata_df (DataFrame): metadata数据框
        output_file (str): 输出h5ad文件路径
        max_bins (int): 最大bin数
    """
    print(f"开始处理Stage，共{len(stage_files)}个细胞文件")
    
    # 收集所有细胞的接触矩阵
    all_matrices = []
    cell_names = []
    cell_metadata = []
    
    # 处理每个细胞文件
    for i, file_info in enumerate(stage_files):
        filename = file_info['filename']
        filepath = file_info['filepath']
        cellname = file_info['cellname']
        
        if i % 100 == 0:
            print(f"  处理进度: {i}/{len(stage_files)} - {filename}")
        
        try:
            # 解析pairs文件
            contacts = parse_pairs_file(filepath)
            
            # 创建稀疏矩阵
            matrix = create_sparse_matrix(contacts, max_bins)
            
            # 保存矩阵和元数据
            all_matrices.append(matrix)
            cell_names.append(cellname)
            
            # 获取细胞的metadata
            metadata_row = metadata_df[metadata_df['Cellname'] == cellname].iloc[0]
            cell_metadata.append(metadata_row.to_dict())
            
        except Exception as e:
            print(f"  处理文件 {filename} 时出错: {e}")
            continue
    
    if not all_matrices:
        print("  没有成功处理任何细胞文件")
        return
    
    print(f"  成功处理 {len(all_matrices)} 个细胞")
    
    # 合并所有矩阵
    # 由于scipy的vstack对于大型稀疏矩阵可能很慢，我们使用另一种方法
    try:
        # 将所有矩阵转换为csr格式并垂直堆叠
        combined_matrix = csr_matrix(np.vstack([m.toarray() for m in all_matrices]))
    except MemoryError:
        print("  内存不足，使用稀疏矩阵堆叠")
        from scipy.sparse import vstack
        combined_matrix = vstack(all_matrices)
    
    # 创建obs数据框
    obs_df = pd.DataFrame(cell_metadata)
    obs_df.index = cell_names
    
    # 创建var数据框 (bin信息)
    var_df = pd.DataFrame(index=[f"bin_{i}" for i in range(max_bins)])
    
    # 创建AnnData对象
    adata = sc.AnnData(X=combined_matrix, obs=obs_df, var=var_df)
    
    # 保存到h5ad文件
    adata.write(output_file)
    print(f"  Stage数据已保存到: {output_file}")


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
        print(f"处理Stage: {stage}")
        
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