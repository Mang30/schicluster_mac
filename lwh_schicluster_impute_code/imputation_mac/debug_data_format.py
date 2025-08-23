#!/usr/bin/env python3
"""
调试数据格式问题，测试单个文件的插补
"""

import os
import sys
import pandas as pd
import numpy as np
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def analyze_data_file(file_path):
    """分析数据文件格式"""
    logging.info(f"📁 分析文件: {file_path}")
    
    # 读取前几行
    with open(file_path, 'r') as f:
        lines = [f.readline().strip() for _ in range(5)]
    
    logging.info("前5行数据:")
    for i, line in enumerate(lines):
        logging.info(f"  {i+1}: {line}")
    
    # 用 pandas 读取
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
        logging.info(f"📊 数据形状: {df.shape}")
        logging.info(f"📊 列数: {df.shape[1]}")
        logging.info(f"📊 前几行:\n{df.head()}")
        return df
    except Exception as e:
        logging.error(f"❌ 读取失败: {e}")
        return None

def test_schicluster_import():
    """测试 schicluster 导入"""
    try:
        sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')
        from schicluster.impute.impute_chromosome import impute_chromosome, MPS_AVAILABLE
        logging.info(f"✅ schicluster 导入成功, MPS 可用: {MPS_AVAILABLE}")
        return True
    except Exception as e:
        logging.error(f"❌ schicluster 导入失败: {e}")
        return False

def test_impute_chromosome():
    """测试 impute_chromosome 函数"""
    if not test_schicluster_import():
        return False
    
    # 测试文件
    test_file = "/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output/E70/GasfE703153/GasfE703153_chr1.txt"
    chrom_size_path = "/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes.txt"
    output_file = "/tmp/test_output.npz"
    
    logging.info(f"🧪 测试插补: {test_file}")
    
    # 分析数据格式
    df = analyze_data_file(test_file)
    if df is None:
        return False
    
    try:
        from schicluster.impute.impute_chromosome import impute_chromosome
        
        # 使用正确的参数调用
        impute_chromosome(
            chrom='chr1',
            resolution=20000,
            output_path=output_file,
            contact_path=test_file,
            chrom_size_path=chrom_size_path,
            use_mps=True,
            chrom1=0,  # 第一列是 bin1 (0-based)
            pos1=0,    # 第一列是 bin1
            chrom2=0,  # 第二列是 bin2 (但这是同一染色体)
            pos2=1     # 第二列是 bin2
        )
        
        logging.info("✅ 插补测试成功!")
        return True
        
    except Exception as e:
        logging.error(f"❌ 插补测试失败: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    logging.info("🚀 开始数据格式调试")
    
    # 测试数据文件
    test_file = "/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output/E70/GasfE703153/GasfE703153_chr1.txt"
    
    if os.path.exists(test_file):
        df = analyze_data_file(test_file)
        test_impute_chromosome()
    else:
        logging.error(f"❌ 测试文件不存在: {test_file}")
