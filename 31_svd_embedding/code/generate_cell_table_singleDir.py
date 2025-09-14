#!/usr/bin/env python3
"""
生成 hicluster embedding 所需的 cell_table.tsv 文件
该文件包含两列: cell_uid 和 file_path，以制表符分隔，无表头
"""

import os
import glob

def generate_cell_table():
    # 输入目录
    input_dir = "/Volumes/SumSung500/CSU/0_HiRES/hires_data_processing/outputs/20length/impute/100K/chunk0"

    # 输出文件路径
    output_file = "/Volumes/SumSung500/CSU/0_HiRES/31_svd_embedding/code/cell_table.tsv"

    # 查找所有 .cool 文件
    cool_files = []
    for file in os.listdir(input_dir):
        if file.endswith('.cool'):
            cool_files.append(file)

    # 排序确保顺序一致
    cool_files.sort()

    print(f"找到 {len(cool_files)} 个 .cool 文件")

    # 生成 cell_table.tsv
    with open(output_file, 'w') as f:
        for cool_file in cool_files:
            # 提取细胞名称（去掉 .cool 后缀）
            cell_name = cool_file.replace('.cool', '')

            # 构造完整路径
            full_path = os.path.join(input_dir, cool_file)

            # 写入文件：cell_name \t file_path
            f.write(f"{cell_name}\t{full_path}\n")

    print(f"生成的 cell_table.tsv 保存在: {output_file}")
    print(f"总共处理了 {len(cool_files)} 个细胞")

if __name__ == "__main__":
    generate_cell_table()