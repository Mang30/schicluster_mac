#!/usr/bin/env python3
"""
专门处理已转换的 bin 格式 Hi-C 数据的插补函数
"""

import os
import sys
import time
import logging
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

# 添加 schicluster 到路径
sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')

def impute_bin_matrix(bin_contact_file, output_path, chrom_size_file, chrom, 
                     resolution=20000, use_mps=True):
    """
    直接处理 bin 格式的接触矩阵进行插补
    
    Parameters:
    -----------
    bin_contact_file : str
        bin1 bin2 count 格式的接触文件
    output_path : str
        输出文件路径
    chrom_size_file : str
        染色体大小文件
    chrom : str
        染色体名称
    resolution : int
        分辨率
    use_mps : bool
        是否使用 MPS 加速
    """
    from schicluster.impute.impute_chromosome import (
        random_walk_auto, gaussian_filter_mps, MPS_AVAILABLE
    )
    from scipy.sparse import diags, save_npz
    from scipy.ndimage import gaussian_filter
    
    logging.info(f"🧬 开始处理 {chrom} (bin格式)")
    
    # 读取染色体大小
    chrom_sizes = pd.read_csv(chrom_size_file, sep='\t', header=None, index_col=0).squeeze()
    n_bins = (chrom_sizes.loc[chrom] // resolution) + 1
    
    # 读取接触数据 (bin1 bin2 count)
    contacts = pd.read_csv(bin_contact_file, sep='\t', header=None, names=['bin1', 'bin2', 'count'])
    
    # 创建稀疏矩阵
    A = csr_matrix((contacts['count'], (contacts['bin1'], contacts['bin2'])), 
                   shape=(n_bins, n_bins))
    A = A + A.T  # 对称化
    
    logging.info(f"   矩阵大小: {A.shape}, 非零元素: {A.nnz}, 稀疏度: {A.nnz/(A.shape[0]*A.shape[1]):.6f}")
    
    # 插补参数
    pad = 1
    std = 1.0
    rp = 0.5
    tol = 0.01
    
    # 移除对角线
    A = A - diags(A.diagonal())
    
    # Gaussian 卷积
    start_time = time.time()
    if pad > 0:
        if use_mps and MPS_AVAILABLE:
            logging.info("   使用 MPS 加速 Gaussian 卷积")
            A_filtered = gaussian_filter_mps(A, std, pad)
            A = csr_matrix(A_filtered)
        else:
            logging.info("   使用 CPU Gaussian 卷积")
            A = gaussian_filter(A.astype(np.float32).toarray(), std, order=0, mode='mirror', truncate=pad)
            A = csr_matrix(A)
    
    conv_time = time.time() - start_time
    logging.info(f"   卷积耗时: {conv_time:.3f} 秒")
    
    # 移除对角线
    A = A - diags(A.diagonal())
    
    # Random Walk with Restart
    start_time = time.time()
    B = A + diags((A.sum(axis=0).A.ravel() == 0).astype(int))
    d = diags(1 / B.sum(axis=0).A.ravel())
    P = d.dot(B).astype(np.float32)
    
    E = random_walk_auto(P, rp, tol, use_mps)
    rwr_time = time.time() - start_time
    logging.info(f"   RWR 耗时: {rwr_time:.3f} 秒")
    
    # 归一化
    start_time = time.time()
    E += E.T
    d = E.sum(axis=0).A.ravel()
    d[d == 0] = 1
    b = diags(1 / np.sqrt(d))
    E = b.dot(E).dot(b)
    norm_time = time.time() - start_time
    logging.info(f"   归一化耗时: {norm_time:.3f} 秒")
    
    # 保存结果
    save_npz(output_path, E)
    
    total_time = conv_time + rwr_time + norm_time
    logging.info(f"   ✅ 完成，总耗时: {total_time:.3f} 秒")
    
    return True

# 测试函数
if __name__ == "__main__":
    # 测试单个文件
    test_file = "/Volumes/SumSung500/CSU/0_HiRES/converted_matrices_by_stage_20k/E70/GasdE701001/GasdE701001_chr1.txt"
    output_file = "/Volumes/SumSung500/CSU/0_HiRES/test_output.npz"
    chrom_size_file = "/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes.txt"
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    try:
        success = impute_bin_matrix(test_file, output_file, chrom_size_file, 'chr1', use_mps=True)
        if success:
            print("🎉 测试成功！")
        else:
            print("❌ 测试失败")
    except Exception as e:
        print(f"❌ 测试出错: {e}")
        import traceback
        traceback.print_exc()
