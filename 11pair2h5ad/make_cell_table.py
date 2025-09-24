#!/usr/bin/env python3
"""
基于染色体数量信息生成分离的细胞表文件

只依赖 chromosome_counts_final.csv 文件和数据目录，根据细胞名和文件名信息
自动构建文件路径，生成分离的细胞表文件。

作者: Claude Code
"""

import os
import pandas as pd
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 硬编码路径配置
CHROMOSOME_COUNTS_FILE = "/Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/chromosome_counts_final.csv"
DATA_ROOT = "/Volumes/SumSung500/CSU/0_HiRES/data/hires/GSE223917_RAW"
OUTPUT_DIR = "/Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad"


def construct_file_path(filename, cellname, data_root):
    """
    根据文件名和细胞名构造文件路径

    Args:
        filename: GSM文件名 (如 GSM6998595)
        cellname: 细胞名 (如 GasaE751001)
        data_root: 数据根目录

    Returns:
        构造的文件路径，如果文件存在的话
    """
    # 方法1: 标准格式 GSM_cellname.pairs.gz
    constructed_path = os.path.join(data_root, f"{filename}_{cellname}.pairs.gz")
    if os.path.exists(constructed_path):
        return constructed_path

    # 方法2: 检查是否有其他命名格式的文件
    # 列出目录中所有包含 filename 或 cellname 的文件
    try:
        files_in_dir = os.listdir(data_root)
        for file in files_in_dir:
            if file.endswith('.pairs.gz'):
                # 检查文件名是否包含 GSM ID 和 cellname
                if filename in file and cellname in file:
                    full_path = os.path.join(data_root, file)
                    return full_path
    except Exception as e:
        logger.debug(f"无法列出目录 {data_root}: {e}")

    return None


def main():
    logger.info("="*60)
    logger.info("基于染色体信息生成分离的细胞表文件")
    logger.info("="*60)

    # 检查输入文件
    if not os.path.exists(CHROMOSOME_COUNTS_FILE):
        logger.error(f"染色体计数文件不存在: {CHROMOSOME_COUNTS_FILE}")
        return

    if not os.path.exists(DATA_ROOT):
        logger.error(f"数据目录不存在: {DATA_ROOT}")
        return

    # 读取染色体计数信息
    logger.info("读取染色体计数信息...")
    chr_counts = pd.read_csv(CHROMOSOME_COUNTS_FILE)
    logger.info(f"加载了 {len(chr_counts)} 个细胞的染色体信息")

    # 统计各类细胞数量
    count_stats = chr_counts['chromosome_count'].value_counts().sort_index()
    logger.info("染色体数量统计:")
    for count, num_cells in count_stats.items():
        gender = "雌性" if count == 20 else "雄性" if count == 21 else "未知"
        logger.info(f"  {count} 条染色体 ({gender}): {num_cells} 个细胞")

    # 构建文件路径并验证文件存在性
    logger.info("构建文件路径...")
    valid_cells = []
    missing_files = []

    for _, row in chr_counts.iterrows():
        filename = row['filename']
        cellname = row['cellname']
        chr_count = row['chromosome_count']

        # 跳过 NaN 值
        if pd.isna(cellname) or pd.isna(filename):
            logger.debug(f"跳过无效记录: filename={filename}, cellname={cellname}")
            continue

        # 构造文件路径
        file_path = construct_file_path(filename, cellname, DATA_ROOT)

        if file_path is not None:
            valid_cells.append({
                'cellname': cellname,
                'file_path': file_path,
                'chromosome_count': chr_count,
                'filename': filename,
                'stage': row.get('stage', ''),
                'celltype': row.get('celltype', '')
            })
        else:
            missing_files.append({
                'filename': filename,
                'cellname': cellname,
                'chr_count': chr_count
            })
            logger.debug(f"无法找到文件: {filename}_{cellname}.pairs.gz")

    logger.info(f"成功找到 {len(valid_cells)} 个细胞的数据文件")
    if missing_files:
        logger.warning(f"无法找到 {len(missing_files)} 个细胞的数据文件")

    # 转换为DataFrame
    merged_df = pd.DataFrame(valid_cells)

    # 分离20条和21条染色体的细胞
    cells_20_chr = merged_df[merged_df['chromosome_count'] == 20]
    cells_21_chr = merged_df[merged_df['chromosome_count'] == 21]

    logger.info(f"20条染色体细胞 (雌性): {len(cells_20_chr)} 个")
    logger.info(f"21条染色体细胞 (雄性): {len(cells_21_chr)} 个")

    # 创建输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 生成20条染色体细胞表
    output_20 = os.path.join(OUTPUT_DIR, '20length_cell_table.tsv')
    cells_20_output = cells_20_chr[['cellname', 'file_path']]
    cells_20_output.to_csv(output_20, sep='\t', header=False, index=False)
    logger.info(f"已生成 20 条染色体细胞表: {output_20}")

    # 生成21条染色体细胞表
    output_21 = os.path.join(OUTPUT_DIR, '21length_cell_table.tsv')
    cells_21_output = cells_21_chr[['cellname', 'file_path']]
    cells_21_output.to_csv(output_21, sep='\t', header=False, index=False)
    logger.info(f"已生成 21 条染色体细胞表: {output_21}")

    # 生成详细信息文件（包含阶段和细胞类型信息）
    if len(merged_df) > 0 and 'stage' in merged_df.columns and 'celltype' in merged_df.columns:
        # 20条染色体详细信息
        output_20_detailed = os.path.join(OUTPUT_DIR, '20length_cell_table_detailed.tsv')
        cells_20_detailed = cells_20_chr[['cellname', 'file_path', 'stage', 'celltype']]
        cells_20_detailed.to_csv(output_20_detailed, sep='\t', header=True, index=False)
        logger.info(f"已生成 20 条染色体详细细胞表: {output_20_detailed}")

        # 21条染色体详细信息
        output_21_detailed = os.path.join(OUTPUT_DIR, '21length_cell_table_detailed.tsv')
        cells_21_detailed = cells_21_chr[['cellname', 'file_path', 'stage', 'celltype']]
        cells_21_detailed.to_csv(output_21_detailed, sep='\t', header=True, index=False)
        logger.info(f"已生成 21 条染色体详细细胞表: {output_21_detailed}")

    # 输出统计信息
    print("\n" + "="*60)
    print("细胞表分离完成")
    print("="*60)
    print(f"📊 总有效细胞数: {len(valid_cells)}")
    print(f"♀️  雌性细胞 (20条染色体): {len(cells_20_chr)}")
    print(f"♂️  雄性细胞 (21条染色体): {len(cells_21_chr)}")
    print(f"📁 输出目录: {OUTPUT_DIR}")
    if missing_files:
        print(f"⚠️  未找到数据文件的细胞: {len(missing_files)}")
    print("="*60)

    # 最终验证：检查生成的文件路径是否都存在
    logger.info("验证生成的文件路径...")
    path_errors = 0
    for _, row in merged_df.iterrows():
        if not os.path.exists(row['file_path']):
            path_errors += 1
            logger.debug(f"文件不存在: {row['file_path']}")

    if path_errors > 0:
        logger.warning(f"⚠️  {path_errors} 个文件路径验证失败")
    else:
        logger.info("✅ 所有生成的文件路径验证通过")

    # 输出使用建议
    print("\n" + "="*50)
    print("📋 使用建议:")
    print("="*50)
    print("1. 雌性细胞处理 (内存需求较低):")
    print(f"   /Users/wuhaoliu/mamba/envs/schicluster/bin/python pairs_to_h5ad_converter.py \\")
    print(f"       --cell_table {output_20} \\")
    print(f"       --output ../output/female_cells_20chr.h5ad \\")
    print(f"       --batch_size 200")
    print()
    print("2. 雄性细胞处理:")
    print(f"   /Users/wuhaoliu/mamba/envs/schicluster/bin/python pairs_to_h5ad_converter.py \\")
    print(f"       --cell_table {output_21} \\")
    print(f"       --output ../output/male_cells_21chr.h5ad \\")
    print(f"       --batch_size 200")
    print("="*50)


if __name__ == '__main__':
    main()