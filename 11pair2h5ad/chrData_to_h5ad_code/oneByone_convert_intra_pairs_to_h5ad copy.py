###
# 用法1：批量处理所有染色体（原有功能）不提供 -c 或 --chromosome 参数即可。
# python convert_intra_pairs_to_h5ad.py \
#    -i /path/to/output/split_pairs \
#    -o /path/to/output/h5ad_files \
#    -s /path/to/your/genome.chrom.sizes \
#    -w 8
# 用法2：单个染色体处理 通过 -c chr1 参数指定目标染色体。
# python convert_intra_pairs_to_h5ad.py \
#    -i /path/to/output/split_pairs \
#    -o /path/to/output/h5ad_files \
#    -s /path/to/your/genome.chrom.sizes \
#    -w 8 \
#    -c chr1
# 用法3：只处理 chrX 的所有细胞文件 通过 -c chrX 参数指定目标染色体。
# python convert_intra_pairs_to_h5ad.py \
#     -i /path/to/output/split_pairs \
#     -o /path/to/output/h5ad_files \
#     -s /path/to/your/genome.chrom.sizes \
#     -w 8 \
#     -c chrX
# ###


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
【代码2 - 修改版】
功能：将同染色体的 pairs.gz 文件转换为 h5ad 格式。
      新增功能：可通过命令行参数指定只处理单个染色体。
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
import argparse # 引入 argparse 用于处理命令行参数

# --- 全局函数 (保持不变) ---

def get_feature_idx_from_coords(b1, b2, total_bins):
    """根据上三角矩阵的坐标(b1, b2)和总bin数，动态计算一维特征索引。"""
    if b1 > b2: b1, b2 = b2, b1
    b1, b2, total_bins = np.int64(b1), np.int64(b2), np.int64(total_bins)
    return total_bins * b1 - b1 * (b1 - 1) // 2 + b2 - b1

def convert_single_chr_pairs(args):
    """
    处理单个同染色体 pairs.gz 文件的核心函数 (为并行化设计)。
    
    Args:
        args (tuple): 包含 (input_file, chr_sizes_map, output_dir, resolution) 的元组。
        
    Returns:
        str: 生成的 h5ad 文件路径，或在出错时返回 None。
    """
    input_file, chr_sizes_map, output_dir, resolution = args
    
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
        n_bins = (chr_length + resolution - 1) // resolution
        n_features = n_bins * (n_bins + 1) // 2

        # 2. 读取 pairs 文件并统计接触
        feature_counts = defaultdict(int)
        with gzip.open(input_file, 'rt') as f_in:
            for line in f_in:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 5: continue
                
                pos1, pos2 = int(parts[2]), int(parts[4])
                
                b1 = pos1 // resolution
                b2 = pos2 // resolution
                
                if b1 < n_bins and b2 < n_bins:
                    feature_idx = get_feature_idx_from_coords(b1, b2, n_bins)
                    feature_counts[feature_idx] += 1
        
        if not feature_counts:
            return None

        # 3. 创建稀疏矩阵
        indices = list(feature_counts.keys())
        data = list(feature_counts.values())
        contact_matrix_sparse = csr_matrix((data, indices, [0, len(data)]),
                                           shape=(1, n_features),
                                           dtype=np.float32)

        # 4. 创建 AnnData 对象
        obs = pd.DataFrame(index=[cell_id])
        rows, cols = np.triu_indices(n_bins)
        var = pd.DataFrame(index=[f'contact_{i}' for i in range(n_features)])
        var['bin1_id'] = rows
        var['bin2_id'] = cols
        
        adata = ad.AnnData(X=contact_matrix_sparse, obs=obs, var=var)
        adata.uns['resolution'] = resolution
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
    # --- MODIFIED: 使用 argparse 定义命令行参数 ---
    parser = argparse.ArgumentParser(description="将同染色体 pairs.gz 文件高效转换为 h5ad 格式。")
    parser.add_argument('-i', '--input_dir', required=True, help="包含 split pairs.gz 文件的输入目录 (代码1的输出)")
    parser.add_argument('-o', '--output_dir', required=True, help="用于存储生成的 h5ad 文件的输出目录")
    parser.add_argument('-s', '--chr_sizes', required=True, help="基因组染色体大小文件路径 (例如 mm10.chrom.sizes)")
    parser.add_argument('-r', '--resolution', type=int, default=100000, help="分辨率 (默认: 100000)")
    parser.add_argument('-w', '--workers', type=int, default=cpu_count() // 2, help="并行工作进程数 (默认: CPU核心数的一半)")
    # --- NEW: 新增的可选参数，用于指定单个染色体 ---
    parser.add_argument('-c', '--chromosome', type=str, default=None, help="(可选) 只处理指定的单个染色体 (例如: chr1, chrX)。如果未提供，则处理所有染色体。")
    
    args = parser.parse_args()

    print("--- 开始执行【代码2 - 修改版】 ---")
    
    # 检查输入
    if not os.path.isdir(args.input_dir):
        print(f"错误: 输入目录不存在: {args.input_dir}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.chr_sizes):
        print(f"错误: 染色体大小文件不存在: {args.chr_sizes}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"输入目录 (split pairs): {args.input_dir}")
    print(f"输入染色体大小文件: {args.chr_sizes}")
    print(f"输出目录 (h5ad): {args.output_dir}")
    print(f"分辨率: {args.resolution}")
    print(f"并行工作进程数: {args.workers}")

    # 加载染色体大小文件
    try:
        chr_sizes_df = pd.read_csv(args.chr_sizes, sep='\t', header=None, names=['chrom', 'size'])
        chr_sizes_map = pd.Series(chr_sizes_df['size'].values, index=chr_sizes_df['chrom']).to_dict()
        print(f"成功加载 {len(chr_sizes_map)} 条染色体的大小信息。")
    except Exception as e:
        print(f"错误: 读取染色体大小文件失败: {e}", file=sys.stderr)
        sys.exit(1)

    # --- MODIFIED: 根据 --chromosome 参数决定查找哪些文件 ---
    target_chr = args.chromosome
    if target_chr:
        print(f"\n模式: 只处理指定的染色体 -> {target_chr}")
        # 确保用户输入的染色体名存在
        if target_chr not in chr_sizes_map:
            print(f"错误: 指定的染色体 '{target_chr}' 不在染色体大小文件中。", file=sys.stderr)
            sys.exit(1)
        glob_pattern = os.path.join(args.input_dir, f"*_{target_chr}.pairs.gz")
    else:
        print("\n模式: 处理所有找到的同染色体文件")
        glob_pattern = os.path.join(args.input_dir, "*_chr*.pairs.gz")

    input_files = glob.glob(glob_pattern)
    
    if not input_files:
        print(f"警告: 在输入目录中未找到任何与模式 '{os.path.basename(glob_pattern)}' 匹配的文件。")
        sys.exit(0)
    print(f"共发现 {len(input_files)} 个文件待处理。")
    
    # 准备并行任务列表
    tasks = [(f, chr_sizes_map, args.output_dir, args.resolution) for f in input_files]

    # 使用多进程池执行任务
    print("\n开始并行处理...")
    with Pool(args.workers) as pool:
        results = list(tqdm(pool.imap_unordered(convert_single_chr_pairs, tasks), total=len(tasks), desc="转换文件为h5ad"))
    
    successful_jobs = sum(1 for r in results if r is not None)
    print(f"\n处理完成！成功转换 {successful_jobs} / {len(tasks)} 个文件为 .h5ad 格式。")
    print("--- 【代码2 - 修改版】执行完毕 ---")

if __name__ == '__main__':
    main()