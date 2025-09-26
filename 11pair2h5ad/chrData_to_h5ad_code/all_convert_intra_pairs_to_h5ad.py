#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
【代码2】
功能：将同染色体的 pairs.gz 文件转换为 h5ad 格式。
输入：代码1生成的 _chrN.pairs.gz 文件。
输出：每个输入文件对应一个 .h5ad 文件，内含100kb分辨率的稀疏特征向量。
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import glob

# --- 全局路径配置 ---
# 请根据您的实际路径修改这些变量
INPUT_SPLIT_DIR = "/path/to/output/split_pairs" # 代码1的输出目录
INPUT_CHR_SIZES = "/path/to/your/genome.chrom.sizes" # 基因组染色体大小文件
OUTPUT_H5AD_DIR = "/path/to/output/h5ad_files"   # 存储生成的h5ad文件的目录
RESOLUTION = 100000
NUM_WORKERS = max(1, cpu_count() // 2)

def get_feature_idx_from_coords(b1, b2, total_bins):
    """根据上三角矩阵的坐标(b1, b2)和总bin数，动态计算一维特征索引。"""
    if b1 > b2: b1, b2 = b2, b1
    b1, b2, total_bins = np.int64(b1), np.int64(b2), np.int64(total_bins)
    return total_bins * b1 - b1 * (b1 - 1) // 2 + b2 - b1

def convert_single_chr_pairs(args):
    """
    处理单个同染色体 pairs.gz 文件的核心函数 (为并行化设计)。
    
    Args:
        args (tuple): 包含 (input_file, chr_sizes_map, output_dir) 的元组。
        
    Returns:
        str: 生成的 h5ad 文件路径，或在出错时返回 None。
    """
    input_file, chr_sizes_map, output_dir = args
    
    try:
        # 从文件名解析 cell_id 和染色体名
        basename = os.path.basename(input_file)
        parts = basename.replace('.pairs.gz', '').rsplit('_', 1)
        if len(parts) != 2:
            print(f"警告: 文件名格式不规范，跳过 {basename}", file=sys.stderr)
            return None
        cell_id, chr_name = parts
        
        if chr_name not in chr_sizes_map:
            print(f"警告: 在染色体大小文件中找不到 {chr_name}，跳过 {basename}", file=sys.stderr)
            return None

        # 1. 计算当前染色体的 bin 数量和总特征数
        chr_length = chr_sizes_map[chr_name]
        n_bins = (chr_length + RESOLUTION - 1) // RESOLUTION
        n_features = n_bins * (n_bins + 1) // 2

        # 2. 读取 pairs 文件并统计接触
        feature_counts = defaultdict(int)
        with gzip.open(input_file, 'rt') as f_in:
            for line in f_in:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 5: continue
                
                pos1, pos2 = int(parts[2]), int(parts[4])
                
                b1 = pos1 // RESOLUTION
                b2 = pos2 // RESOLUTION
                
                # 确保坐标在有效范围内
                if b1 < n_bins and b2 < n_bins:
                    feature_idx = get_feature_idx_from_coords(b1, b2, n_bins)
                    feature_counts[feature_idx] += 1
        
        if not feature_counts:
            # print(f"信息: 文件 {basename} 中无有效接触数据，跳过。")
            return None

        # 3. 创建稀疏矩阵
        indices = list(feature_counts.keys())
        data = list(feature_counts.values())
        contact_matrix_sparse = csr_matrix((data, indices, [0, len(data)]),
                                           shape=(1, n_features),
                                           dtype=np.float32)

        # 4. 创建 AnnData 对象
        obs = pd.DataFrame(index=[cell_id])
        
        # 创建包含详细坐标信息的 var DataFrame
        rows, cols = np.triu_indices(n_bins)
        var = pd.DataFrame(index=[f'contact_{i}' for i in range(n_features)])
        var['bin1_id'] = rows
        var['bin2_id'] = cols
        
        adata = ad.AnnData(X=contact_matrix_sparse, obs=obs, var=var)
        adata.uns['resolution'] = RESOLUTION
        adata.uns['chromosome'] = chr_name

        # 5. 保存 h5ad 文件
        output_filename = f"{cell_id}_{chr_name}.h5ad"
        output_path = os.path.join(output_dir, output_filename)
        adata.write_h5ad(output_path)
        
        return output_path

    except Exception as e:
        print(f"错误: 处理文件 {input_file} 时发生错误: {e}", file=sys.stderr)
        return None

def main():
    """主执行函数"""
    print("--- 开始执行【代码2】：将同染色体 pairs.gz 转换为 h5ad ---")
    
    # 检查输入
    if not os.path.isdir(INPUT_SPLIT_DIR):
        print(f"错误: 输入目录不存在: {INPUT_SPLIT_DIR}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(INPUT_CHR_SIZES):
        print(f"错误: 染色体大小文件不存在: {INPUT_CHR_SIZES}", file=sys.stderr)
        sys.exit(1)

    # 创建输出目录
    os.makedirs(OUTPUT_H5AD_DIR, exist_ok=True)
    print(f"输入目录 (split pairs): {INPUT_SPLIT_DIR}")
    print(f"输入染色体大小文件: {INPUT_CHR_SIZES}")
    print(f"输出目录 (h5ad): {OUTPUT_H5AD_DIR}")
    print(f"分辨率: {RESOLUTION}")
    print(f"并行工作进程数: {NUM_WORKERS}")

    # 加载染色体大小文件
    try:
        chr_sizes_df = pd.read_csv(INPUT_CHR_SIZES, sep='\t', header=None, names=['chrom', 'size'])
        chr_sizes_map = pd.Series(chr_sizes_df['size'].values, index=chr_sizes_df['chrom']).to_dict()
        print(f"成功加载 {len(chr_sizes_map)} 条染色体的大小信息。")
    except Exception as e:
        print(f"错误: 读取染色体大小文件失败: {e}", file=sys.stderr)
        sys.exit(1)

    # 查找所有待处理的同染色体文件
    input_files = glob.glob(os.path.join(INPUT_SPLIT_DIR, "*_chr*.pairs.gz"))
    if not input_files:
        print("警告: 在输入目录中未找到任何 `_chr*.pairs.gz` 文件。")
        sys.exit(0)
    print(f"共发现 {len(input_files)} 个同染色体文件待处理。")
    
    # 准备并行任务列表
    tasks = [(f, chr_sizes_map, OUTPUT_H5AD_DIR) for f in input_files]

    # 使用多进程池执行任务
    print("\n开始并行处理...")
    with Pool(NUM_WORKERS) as pool:
        results = list(tqdm(pool.imap_unordered(convert_single_chr_pairs, tasks), total=len(tasks), desc="转换文件为h5ad"))
    
    successful_jobs = sum(1 for r in results if r is not None)
    print(f"\n处理完成！成功转换 {successful_jobs} / {len(tasks)} 个文件为 .h5ad 格式。")
    print("--- 【代码2】执行完毕 ---")

if __name__ == '__main__':
    main()