#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
细胞类型注释匹配逻辑实现
"""

import pandas as pd
import os
import re
from collections import defaultdict


def extract_cellname_from_filename(filename):
    """
    从scHi-C文件名中提取Cellname
    
    Args:
        filename (str): scHi-C文件名，格式如 GSM6998595_GasaE751001.pairs.gz
        
    Returns:
        str: 提取的Cellname，如 GasaE751001
    """
    # 匹配 _ 之后 .pairs.gz 之前的字符串
    match = re.search(r'_([^_]+)\.pairs\.gz$', filename)
    if match:
        return match.group(1)
    return None


def match_cell_metadata(data_dir, metadata_file):
    """
    匹配scHi-C文件和metadata中的细胞信息
    
    Args:
        data_dir (str): scHi-C数据目录路径
        metadata_file (str): metadata文件路径
        
    Returns:
        dict: 按Stage分组的文件信息
    """
    # 读取metadata
    metadata_df = pd.read_excel(metadata_file)
    
    # 创建Cellname到metadata的映射
    cell_metadata = {}
    for _, row in metadata_df.iterrows():
        cell_metadata[row['Cellname']] = row.to_dict()
    
    # 获取所有scHi-C文件
    pairs_files = [f for f in os.listdir(data_dir) if f.endswith('.pairs.gz')]
    
    # 按Stage分组存储文件信息
    stage_files = defaultdict(list)
    
    # 匹配文件和metadata
    matched_count = 0
    unmatched_files = []
    
    for filename in pairs_files:
        cellname = extract_cellname_from_filename(filename)
        if cellname and cellname in cell_metadata:
            # 获取细胞的metadata信息
            cell_info = cell_metadata[cellname]
            stage = cell_info['Stage']
            
            # 存储文件路径和细胞信息
            file_info = {
                'filename': filename,
                'filepath': os.path.join(data_dir, filename),
                'cellname': cellname,
                'stage': stage,
                'celltype': cell_info['Celltype']
            }
            
            stage_files[stage].append(file_info)
            matched_count += 1
        else:
            unmatched_files.append(filename)
    
    print(f"匹配成功的文件数: {matched_count}")
    print(f"未匹配的文件数: {len(unmatched_files)}")
    
    if unmatched_files:
        print("前10个未匹配的文件:")
        for f in unmatched_files[:10]:
            print(f"  {f}")
    
    # 显示各Stage的文件数量
    print("\n各Stage的文件数量:")
    for stage, files in stage_files.items():
        print(f"  {stage}: {len(files)} 个文件")
    
    return stage_files


def main():
    """主函数"""
    data_dir = "/home/duxuyan/Projects/0_HiRES/data/GSE223917_RAW"
    metadata_file = "/home/duxuyan/Projects/0_HiRES/data/GSE223917_HiRES_emb_metadata.xlsx"
    
    # 执行匹配
    stage_files = match_cell_metadata(data_dir, metadata_file)
    
    # 保存匹配结果到文件
    output_file = "/home/duxuyan/Projects/0_HiRES/code/claude_new/stage_files_mapping.csv"
    
    # 展平数据用于保存
    all_files_info = []
    for stage, files in stage_files.items():
        for file_info in files:
            all_files_info.append(file_info)
    
    # 转换为DataFrame并保存
    df = pd.DataFrame(all_files_info)
    df.to_csv(output_file, index=False)
    print(f"\n匹配结果已保存到: {output_file}")


if __name__ == "__main__":
    main()