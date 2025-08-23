#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
继续处理剩余Stage的h5ad文件生成脚本
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys


def process_single_stage(stage_name):
    """处理单个Stage"""
    print(f"开始处理Stage: {stage_name}")
    
    # 添加代码路径
    sys.path.append('/home/duxuyan/Projects/0_HiRES/code/claude_new')
    from h5ad_generation_fixed import process_stage_to_h5ad
    
    # 读取stage文件映射
    stage_files_mapping = "/home/duxuyan/Projects/0_HiRES/code/claude_new/stage_files_mapping.csv"
    metadata_file = "/home/duxuyan/Projects/0_HiRES/data/GSE223917_HiRES_emb_metadata.xlsx"
    
    # 读取数据
    stage_df = pd.read_csv(stage_files_mapping)
    metadata_df = pd.read_excel(metadata_file)
    
    # 获取该Stage的所有文件
    stage_files = stage_df[stage_df['stage'] == stage_name].to_dict('records')
    
    print(f"Stage {stage_name} 包含 {len(stage_files)} 个细胞文件")
    
    # 定义输出文件路径
    output_dir = "/home/duxuyan/Projects/0_HiRES/output"
    output_file = os.path.join(output_dir, f"stage_{stage_name}.h5ad")
    
    # 处理该Stage的数据
    process_stage_to_h5ad(stage_files, metadata_df, output_file, max_bins=2000)


def main():
    """主函数"""
    # 处理剩余的Stages
    remaining_stages = ['EX15', 'E95', 'EX05', 'E85', 'E75']
    
    for stage in remaining_stages:
        try:
            # 检查是否已经处理过
            output_file = f"/home/duxuyan/Projects/0_HiRES/output/stage_{stage}.h5ad"
            if os.path.exists(output_file):
                print(f"Stage {stage} 已经处理过，跳过")
                continue
                
            process_single_stage(stage)
        except Exception as e:
            print(f"处理Stage {stage} 时出错: {e}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()