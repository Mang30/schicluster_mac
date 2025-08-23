#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
测试Stage处理脚本 - 用于验证修复后的h5ad生成逻辑
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
from scipy.sparse import csr_matrix


def test_feature_vector_creation():
    """测试特征向量创建功能"""
    print("测试特征向量创建功能...")
    
    # 创建测试数据
    contacts = {
        (0, 0): 5,
        (0, 1): 3,
        (1, 1): 2,
        (2, 3): 4,
        (3, 3): 1
    }
    
    # 导入修复后的函数
    sys.path.append('/home/duxuyan/Projects/0_HiRES/code/claude_new')
    from h5ad_generation_fixed import create_feature_vector
    
    # 创建特征向量
    vector = create_feature_vector(contacts, max_bins=5)
    
    print(f"特征向量形状: {vector.shape}")
    print(f"非零元素数量: {vector.nnz}")
    print(f"非零元素: {vector.data}")
    print(f"列索引: {vector.indices}")
    
    # 验证结果
    expected_shape = (1, 15)  # 5*6/2 = 15个上三角元素
    assert vector.shape == expected_shape, f"形状不匹配: {vector.shape} != {expected_shape}"
    assert vector.nnz > 0, "没有非零元素"
    
    print("特征向量创建测试通过!")


def test_small_stage_processing():
    """测试小规模Stage处理"""
    print("\n测试小规模Stage处理...")
    
    # 读取stage文件映射
    stage_files_mapping = "/home/duxuyan/Projects/0_HiRES/code/claude_new/stage_files_mapping.csv"
    metadata_file = "/home/duxuyan/Projects/0_HiRES/data/GSE223917_HiRES_emb_metadata.xlsx"
    
    if not os.path.exists(stage_files_mapping):
        print(f"文件不存在: {stage_files_mapping}")
        return
        
    # 读取数据
    stage_df = pd.read_csv(stage_files_mapping)
    metadata_df = pd.read_excel(metadata_file)
    
    # 选择最小的Stage进行测试
    stage_counts = stage_df['stage'].value_counts()
    smallest_stage = stage_counts.idxmin()
    print(f"选择最小Stage进行测试: {smallest_stage} ({stage_counts.min()} 个细胞)")
    
    # 获取该Stage的前5个文件进行测试
    stage_files = stage_df[stage_df['stage'] == smallest_stage].to_dict('records')[:5]
    
    print(f"测试处理 {len(stage_files)} 个细胞文件")
    
    # 导入处理函数
    sys.path.append('/home/duxuyan/Projects/0_HiRES/code/claude_new')
    from h5ad_generation_fixed import process_stage_to_h5ad
    
    # 创建输出目录
    output_dir = "/home/duxuyan/Projects/0_HiRES/output/test"
    os.makedirs(output_dir, exist_ok=True)
    
    # 定义输出文件路径
    output_file = os.path.join(output_dir, f"test_stage_{smallest_stage}.h5ad")
    
    # 处理该Stage的数据
    process_stage_to_h5ad(stage_files, metadata_df, output_file, max_bins=1000)
    
    # 验证生成的文件
    if os.path.exists(output_file):
        print(f"测试文件生成成功: {output_file}")
        
        # 读取并验证文件
        adata = sc.read_h5ad(output_file)
        print(f"生成的AnnData对象形状: {adata.shape}")
        print(f"obs行数: {adata.obs.shape[0]}")
        print(f"var行数: {adata.var.shape[0]}")
        
        # 验证维度一致性
        assert adata.shape[0] == len(stage_files), f"行数不匹配: {adata.shape[0]} != {len(stage_files)}"
        print("文件验证通过!")
    else:
        print("测试文件生成失败")


def main():
    """主函数"""
    print("开始测试修复后的h5ad生成逻辑...")
    
    try:
        # 运行测试
        test_feature_vector_creation()
        test_small_stage_processing()
        
        print("\n所有测试完成!")
    except Exception as e:
        print(f"测试过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()