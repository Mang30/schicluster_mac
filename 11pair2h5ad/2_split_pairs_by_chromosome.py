#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
【代码1】
功能：将单个细胞的全基因组 pairs.gz 文件按染色体拆分。
输入：cell_table.tsv (含 cell_id 和 pairs.gz 文件路径)。
输出：每个细胞对应一组拆分后的文件，包括
      - {cell_id}_chrN.pairs.gz (N为染色体号)
      - {cell_id}_interchrom.pairs.gz (跨染色体)
"""
import os
import sys
import gzip
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# --- 全局路径配置 ---
# 请根据您的实际路径修改这些变量
INPUT_CELL_TABLE = "/home/duxuyan/Projects/schicluster_mac/11pair2h5ad/cell_table.tsv"
OUTPUT_SPLIT_DIR = "/home/duxuyan/Projects/schicluster_mac/11pair2h5ad/output/split_pairs" # 存储拆分后文件的目录
NUM_WORKERS = max(1, cpu_count() // 2)  # 使用一半的CPU核心进行并行处理，可按需调整

def split_single_pairs_file(args):
    """
    处理单个 pairs.gz 文件的核心函数 (为并行化设计)。
    
    Args:
        args (tuple): 包含 (cell_id, input_pairs_path, output_dir) 的元组。
    
    Returns:
        str: 处理完成的输入文件路径，或在出错时返回 None。
    """
    cell_id, input_pairs_path, output_dir = args
    
    if not os.path.exists(input_pairs_path):
        print(f"警告: 文件不存在，跳过 {input_pairs_path}", file=sys.stderr)
        return None

    output_files = {}  # 存储输出文件的句柄，避免频繁开关文件
    
    try:
        with gzip.open(input_pairs_path, 'rt') as f_in:
            for line in f_in:
                # 跳过注释行
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                # .pairs 格式通常至少有 chr1, pos1, chr2, pos2 四列
                if len(parts) < 4:
                    continue
                    
                chr1, chr2 = parts[1], parts[3]
                
                # 判断是同染色体还是跨染色体
                if chr1 == chr2:
                    output_key = chr1
                else:
                    output_key = "interchrom"
                
                # 如果是第一次遇到这个key，则创建并打开对应的输出文件
                if output_key not in output_files:
                    output_filename = f"{cell_id}_{output_key}.pairs.gz"
                    output_path = os.path.join(output_dir, output_filename)
                    output_files[output_key] = gzip.open(output_path, 'wt')
                
                # 写入数据
                output_files[output_key].write(line)
        
        return input_pairs_path

    except Exception as e:
        print(f"错误: 处理文件 {input_pairs_path} 时发生错误: {e}", file=sys.stderr)
        return None
    
    finally:
        # 确保所有打开的文件都被关闭
        for f_handle in output_files.values():
            f_handle.close()

def main():
    """主执行函数"""
    print("--- 开始执行【代码1】：按染色体拆分 pairs.gz 文件 ---")
    
    # 检查输入文件
    if not os.path.exists(INPUT_CELL_TABLE):
        print(f"错误: 输入的 cell_table 文件不存在: {INPUT_CELL_TABLE}", file=sys.stderr)
        sys.exit(1)
        
    # 创建输出目录
    os.makedirs(OUTPUT_SPLIT_DIR, exist_ok=True)
    print(f"输入 cell_table: {INPUT_CELL_TABLE}")
    print(f"输出目录: {OUTPUT_SPLIT_DIR}")
    print(f"并行工作进程数: {NUM_WORKERS}")

    # 加载细胞信息表
    try:
        cell_df = pd.read_csv(INPUT_CELL_TABLE, sep='\t', header=None, names=['cell_id', 'file_path'])
        print(f"共发现 {len(cell_df)} 个细胞待处理。")
    except Exception as e:
        print(f"错误: 读取 cell_table 文件失败: {e}", file=sys.stderr)
        sys.exit(1)

    # 准备并行任务列表
    tasks = []
    for _, row in cell_df.iterrows():
        tasks.append((row['cell_id'], row['file_path'], OUTPUT_SPLIT_DIR))

    # 使用多进程池执行任务
    print("\n开始并行处理...")
    with Pool(NUM_WORKERS) as pool:
        # 使用 tqdm 显示进度条
        results = list(tqdm(pool.imap_unordered(split_single_pairs_file, tasks), total=len(tasks), desc="拆分文件"))

    successful_jobs = sum(1 for r in results if r is not None)
    print(f"\n处理完成！成功拆分 {successful_jobs} / {len(tasks)} 个细胞的 pairs 文件。")
    print("--- 【代码1】执行完毕 ---")

if __name__ == '__main__':
    main()