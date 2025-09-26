#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
【代码2 - 重构版】
功能：为指定的单个染色体，聚合所有细胞的 pairs.gz 数据，并生成一个
      包含所有细胞信息的 h5ad 文件。
用法：为每条需要处理的染色体分别运行一次本脚本。
"""

import os
import sys
import gzip
import pandas as pd
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix, vstack
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import glob
import argparse

# --- 全局函数 (保持不变) ---

def get_feature_idx_from_coords(b1, b2, total_bins):
    """根据上三角矩阵的坐标(b1, b2)和总bin数，动态计算一维特征索引。"""
    if b1 > b2: b1, b2 = b2, b1
    b1, b2, total_bins = np.int64(b1), np.int64(b2), np.int64(total_bins)
    return total_bins * b1 - b1 * (b1 - 1) // 2 + b2 - b1

def process_one_cell_for_chromosome(args):
    """
    处理单个细胞关于特定染色体的 pairs.gz 文件。
    不再创建AnnData对象，只返回稀疏向量和细胞ID，用于后续聚合。
    """
    input_file, chr_name, chr_sizes_map, resolution = args
    
    try:
        basename = os.path.basename(input_file)
        cell_id = basename.replace(f'_{chr_name}.pairs.gz', '')
        
        chr_length = chr_sizes_map[chr_name]
        n_bins = (chr_length + resolution - 1) // resolution
        n_features = n_bins * (n_bins + 1) // 2

        feature_counts = defaultdict(int)
        with gzip.open(input_file, 'rt') as f_in:
            for line in f_in:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 5: continue
                
                pos1, pos2 = int(parts[2]), int(parts[4])
                b1, b2 = pos1 // resolution, pos2 // resolution
                
                if b1 < n_bins and b2 < n_bins:
                    feature_idx = get_feature_idx_from_coords(b1, b2, n_bins)
                    feature_counts[feature_idx] += 1
        
        if not feature_counts:
            return None

        indices = list(feature_counts.keys())
        data = list(feature_counts.values())
        contact_matrix_sparse = csr_matrix((data, indices, [0, len(data)]),
                                           shape=(1, n_features),
                                           dtype=np.float32)
        
        return contact_matrix_sparse, cell_id

    except Exception as e:
        print(f"错误: 处理文件 {input_file} 时发生错误: {e}", file=sys.stderr)
        return None

def create_and_save_anndata(X_sparse, cell_ids, chr_name, n_bins, resolution, output_path):
    """
    根据聚合后的数据创建并保存最终的 AnnData 对象。
    """
    logger = logging.getLogger(__name__)
    logger.info(f"为染色体 {chr_name} 创建最终的 AnnData 对象...")

    obs = pd.DataFrame(index=cell_ids)
    
    logger.info("创建特征(var)注释...")
    n_features = n_bins * (n_bins + 1) // 2
    rows, cols = np.triu_indices(n_bins)
    var = pd.DataFrame(index=[f'contact_{i}' for i in range(n_features)])
    var['bin1_id'] = rows
    var['bin2_id'] = cols
    
    adata = ad.AnnData(X=X_sparse, obs=obs, var=var)
    adata.uns['resolution'] = resolution
    adata.uns['chromosome'] = chr_name
    
    logger.info(f"保存聚合后的 h5ad 文件到: {output_path}")
    adata.write_h5ad(output_path)
    
    logger.info("-" * 30 + " 聚合完成! " + "-" * 30)
    logger.info(f"最终数据形状 (细胞数 x 特征数): {adata.shape}")
    logger.info(f"稀疏度: {adata.X.nnz / (adata.shape[0] * adata.shape[1]) * 100:.4f}%")

def main():
    """主执行函数"""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    parser = argparse.ArgumentParser(description="为指定的单个染色体，聚合所有细胞的 pairs 数据并生成一个 h5ad 文件。")
    parser.add_argument('-i', '--input_dir', required=True, help="包含 split pairs.gz 文件的输入目录 (代码1的输出)")
    parser.add_argument('-o', '--output_file', required=True, help="输出的 h5ad 文件名 (例如: /path/to/chr1.h5ad)")
    parser.add_argument('-s', '--chr_sizes', required=True, help="基因组染色体大小文件路径")
    parser.add_argument('-c', '--chromosome', required=True, type=str, help="必须指定要处理的单个染色体 (例如: chr1, chrX)")
    parser.add_argument('-r', '--resolution', type=int, default=100000, help="分辨率 (默认: 100000)")
    parser.add_argument('-w', '--workers', type=int, default=cpu_count() // 2, help="并行工作进程数")
    
    args = parser.parse_args()
    
    logger = logging.getLogger(__name__)
    logger.info(f"--- 开始执行【代码2 - 重构版】，目标染色体: {args.chromosome} ---")

    # 检查输入
    if not os.path.isdir(args.input_dir):
        logger.error(f"错误: 输入目录不存在: {args.input_dir}")
        sys.exit(1)
    if not os.path.exists(args.chr_sizes):
        logger.error(f"错误: 染色体大小文件不存在: {args.chr_sizes}")
        sys.exit(1)
    
    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # 加载染色体大小
    chr_sizes_df = pd.read_csv(args.chr_sizes, sep='\t', header=None, names=['chrom', 'size'])
    chr_sizes_map = pd.Series(chr_sizes_df['size'].values, index=chr_sizes_df['chrom']).to_dict()
    
    target_chr = args.chromosome
    if target_chr not in chr_sizes_map:
        logger.error(f"错误: 指定的染色体 '{target_chr}' 不在染色体大小文件中。")
        sys.exit(1)

    # 查找所有与目标染色体相关的细胞文件
    glob_pattern = os.path.join(args.input_dir, f"*_{target_chr}.pairs.gz")
    input_files = sorted(glob.glob(glob_pattern)) # 排序以保证一定程度的确定性
    
    if not input_files:
        logger.warning(f"警告: 在输入目录中未找到任何与模式 '{os.path.basename(glob_pattern)}' 匹配的文件。")
        sys.exit(0)
    logger.info(f"共发现 {len(input_files)} 个关于 {target_chr} 的细胞文件待处理。")
    
    # 准备并行任务
    tasks = [(f, target_chr, chr_sizes_map, args.resolution) for f in input_files]

    # 并行处理所有细胞文件
    logger.info("开始并行处理所有细胞文件...")
    results = []
    with Pool(args.workers) as pool:
        for res in tqdm(pool.imap_unordered(process_one_cell_for_chromosome, tasks), total=len(tasks), desc=f"处理 {target_chr} 数据"):
            if res is not None:
                results.append(res)
    
    if not results:
        logger.error("错误: 未能成功处理任何细胞数据。")
        sys.exit(1)
    
    logger.info(f"成功处理了 {len(results)} 个细胞的数据，现在开始聚合...")

    # 按细胞ID（字符串）排序，确保每次运行的 h5ad 行顺序都一致
    results.sort(key=lambda x: x[1])

    # 拆分结果为稀疏矩阵列表和细胞ID列表
    cell_matrices_sparse = [res[0] for res in results]
    cell_ids = [res[1] for res in results]

    # 聚合所有稀疏矩阵
    X_aggregated_sparse = vstack(cell_matrices_sparse, format='csr')

    # 计算当前染色体的bin数量
    chr_length = chr_sizes_map[target_chr]
    n_bins = (chr_length + args.resolution - 1) // args.resolution

    # 创建并保存最终的 AnnData 对象
    create_and_save_anndata(X_aggregated_sparse, cell_ids, target_chr, n_bins, args.resolution, args.output_file)
    
    logger.info(f"--- 【代码2 - 重构版】为染色体 {args.chromosome} 的处理全部完成 ---")

if __name__ == '__main__':
    main()
    