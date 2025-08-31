#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path
import logging

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def read_stage_info(csv_file_path):
    df = pd.read_csv(csv_file_path)
    stage_mapping = {}
    for _, row in df.iterrows():
        cellname = row['cellname']
        stage = row['stage']
        stage_mapping[cellname] = stage
    return stage_mapping

def create_stage_directories(output_base_dir, stages):
    stage_dirs = {}
    for stage in stages:
        stage_dir = Path(output_base_dir) / stage
        stage_dir.mkdir(parents=True, exist_ok=True)
        stage_dirs[stage] = stage_dir
        logging.info(f"创建目录: {stage_dir}")
    return stage_dirs

def create_symlinks_by_stage(cool_files_dir, stage_mapping, stage_dirs):
    cool_files_dir = Path(cool_files_dir)
    created_links = 0
    missing_files = []
    
    for cool_file in cool_files_dir.glob("*.cool"):
        cellname = cool_file.stem
        
        if cellname in stage_mapping:
            stage = stage_mapping[cellname]
            if stage in stage_dirs:
                target_link = stage_dirs[stage] / cool_file.name
                
                if target_link.exists():
                    target_link.unlink()
                
                try:
                    os.symlink(cool_file.absolute(), target_link)
                    created_links += 1
                    logging.info(f"创建软链接: {cellname} -> {stage}")
                except Exception as e:
                    logging.error(f"创建软链接失败 {cellname}: {e}")
            else:
                logging.warning(f"未找到stage目录: {stage}")
        else:
            missing_files.append(cellname)
            logging.warning(f"未在CSV中找到细胞: {cellname}")
    
    if missing_files:
        logging.info(f"共有 {len(missing_files)} 个文件未在CSV中找到对应的stage信息")
    
    return created_links, missing_files

def organize_cool_files(cool_files_dir, csv_file_path, output_base_dir):
    logging.info(f"开始处理cool文件组织...")
    logging.info(f"Cool文件目录: {cool_files_dir}")
    logging.info(f"CSV文件路径: {csv_file_path}")
    logging.info(f"输出目录: {output_base_dir}")
    
    stage_mapping = read_stage_info(csv_file_path)
    logging.info(f"从CSV文件读取到 {len(stage_mapping)} 个细胞的stage信息")
    
    unique_stages = set(stage_mapping.values())
    logging.info(f"发现 {len(unique_stages)} 个不同的stage: {sorted(unique_stages)}")
    
    stage_dirs = create_stage_directories(output_base_dir, unique_stages)
    
    created_links, missing_files = create_symlinks_by_stage(
        cool_files_dir, stage_mapping, stage_dirs
    )
    
    logging.info(f"处理完成!")
    logging.info(f"成功创建 {created_links} 个软链接")
    
    return created_links, missing_files

if __name__ == "__main__":
    setup_logging()
    
    # 输入和输出路径配置
    cool_files_20length_dir = "outputs/20length/impute/100K/chunk0"
    cool_files_21length_dir = "outputs/21length/impute/100K/chunk0"
    csv_20length_file = "../20length_chromosome_final_with_path.csv"
    csv_21length_file = "../21length_chromosome_final_with_path.csv"
    output_20length_dir = "outputs/20length/stage_organized"
    output_21length_dir = "outputs/21length/stage_organized"
    
    # 处理20length数据（用于测试）
    logging.info("=" * 50)
    logging.info("处理20length数据")
    logging.info("=" * 50)
    
    created_links_20, missing_files_20 = organize_cool_files(
        cool_files_20length_dir,
        csv_20length_file,
        output_20length_dir
    )
    
    # 可以选择性处理21length数据（当数据准备好时）
    process_21length = False  # 设置为True来处理21length数据
    
    if process_21length and Path(cool_files_21length_dir).exists():
        logging.info("\n" + "=" * 50)
        logging.info("处理21length数据")
        logging.info("=" * 50)
        
        created_links_21, missing_files_21 = organize_cool_files(
            cool_files_21length_dir,
            csv_21length_file,
            output_21length_dir
        )
    else:
        logging.info("\n跳过21length数据处理（数据尚未准备好或未启用）")