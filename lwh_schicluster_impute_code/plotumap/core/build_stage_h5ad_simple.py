#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
简化版本：为一个 stage 构建 h5ad 文件
功能：
 - 收集 stage 下所有单细胞目录
 - 每个细胞提取所有染色体的上三角特征并拼接
 - 构建特征矩阵：行=细胞，列=上三角特征
 - 保存为 h5ad 文件

用法：
  python build_stage_h5ad_simple.py \
    --stage-dir /path/to/E70 \
    --obs-xlsx /path/to/metadata.xlsx \
    --output /path/to/E70.h5ad \
    --max-cells 100
"""

import argparse
import logging
import re
from pathlib import Path
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import sparse
import anndata as ad

logger = logging.getLogger('build_stage_h5ad')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def extract_chr_label(filename: str) -> str:
    """从文件名提取染色体标签"""
    m = re.search(r"(chr(?:[0-9]+|X))", filename)
    return m.group(1) if m else None


def sort_chromosomes(chr_list: List[str]) -> List[str]:
    """对染色体进行排序：1,2,...,19,X"""
    def chr_key(label: str):
        if label == 'chrX':
            return 1000
        try:
            return int(label[3:])
        except:
            return 999
    return sorted(chr_list, key=chr_key)


def discover_chromosome_structure(stage_dir: Path) -> Tuple[List[str], Dict[str, int]]:
    """
    扫描第一个细胞目录，确定染色体结构
    返回：(排序的染色体列表, 每个染色体的bin数)
    """
    cell_dirs = [d for d in stage_dir.iterdir() if d.is_dir() and any(d.glob('*.npz'))]
    if not cell_dirs:
        raise ValueError(f"在 {stage_dir} 中未找到包含 .npz 文件的细胞目录")
    
    # 使用第一个细胞确定结构
    first_cell = cell_dirs[0]
    chr_info = {}
    
    for npz_file in first_cell.glob('*.npz'):
        chr_label = extract_chr_label(npz_file.name)
        if chr_label:
            arr = np.load(npz_file, allow_pickle=True)
            if 'shape' in arr.files:
                shape = tuple(arr['shape'])
                n_bins = int(shape[0])
                chr_info[chr_label] = n_bins
    
    # 排序染色体
    chr_list = sort_chromosomes(list(chr_info.keys()))
    
    logger.info(f"发现 {len(chr_list)} 条染色体: {chr_list}")
    total_features = sum(n * (n + 1) // 2 for n in chr_info.values())
    logger.info(f"总特征数: {total_features:,}")
    
    return chr_list, chr_info


def extract_cell_triu_features(cell_dir: Path, chr_list: List[str], chr_info: Dict[str, int]) -> np.ndarray:
    """
    提取单个细胞所有染色体的上三角特征
    """
    all_features = []
    
    for chr_label in chr_list:
        # 找到对应的 npz 文件
        npz_file = None
        for f in cell_dir.glob('*.npz'):
            if extract_chr_label(f.name) == chr_label:
                npz_file = f
                break
        
        if npz_file is None:
            raise FileNotFoundError(f"细胞 {cell_dir.name} 缺少染色体 {chr_label}")
        
        # 加载并重建稀疏矩阵
        arr = np.load(npz_file, allow_pickle=True)
        csr = sparse.csr_matrix(
            (arr['data'], arr['indices'], arr['indptr']), 
            shape=tuple(arr['shape'])
        )
        
        # 提取上三角
        n = csr.shape[0]
        triu = sparse.triu(csr, k=0).tocoo()
        
        # 转换为线性索引的dense向量
        n_tri = n * (n + 1) // 2
        tri_vec = np.zeros(n_tri, dtype=np.float32)
        
        if len(triu.data) > 0:
            # 计算上三角的线性索引
            row_offsets = (np.arange(n) * n) - (np.arange(n) * (np.arange(n) - 1) // 2)
            linear_idx = row_offsets[triu.row] + (triu.col - triu.row)
            tri_vec[linear_idx] = triu.data.astype(np.float32)
        
        all_features.append(tri_vec)
    
    return np.concatenate(all_features)


def load_metadata(obs_xlsx: Path) -> pd.DataFrame:
    """加载细胞metadata"""
    df = pd.read_excel(obs_xlsx)
    
    # 尝试找到细胞ID列
    id_cols = ['cell', 'cell_id', 'cell_name', 'Cell', 'CellID', 'cellID', 'cellName']
    for col in id_cols:
        if col in df.columns:
            df = df.set_index(col)
            logger.info(f"使用 {col} 作为细胞ID索引")
            break
    else:
        logger.warning("未找到明确的细胞ID列，使用默认索引")
    
    return df


def build_stage_h5ad(stage_dir: Path, obs_xlsx: Path, output_path: Path, max_cells: int = None):
    """
    为一个 stage 构建 h5ad 文件
    """
    logger.info(f"开始处理 stage: {stage_dir}")
    
    # 1. 发现染色体结构
    chr_list, chr_info = discover_chromosome_structure(stage_dir)
    
    # 2. 收集所有细胞目录
    cell_dirs = [d for d in stage_dir.iterdir() if d.is_dir() and any(d.glob('*.npz'))]
    cell_dirs = sorted(cell_dirs, key=lambda x: x.name)
    
    if max_cells:
        cell_dirs = cell_dirs[:max_cells]
    
    logger.info(f"将处理 {len(cell_dirs)} 个细胞")
    
    # 3. 加载metadata
    obs_df = load_metadata(obs_xlsx)
    
    # 4. 提取所有细胞的特征
    X_list = []
    obs_list = []
    failed_cells = []
    
    for cell_dir in tqdm(cell_dirs, desc="提取细胞特征"):
        cell_id = cell_dir.name
        try:
            features = extract_cell_triu_features(cell_dir, chr_list, chr_info)
            X_list.append(features)
            
            # 获取metadata
            if cell_id in obs_df.index:
                cell_obs = obs_df.loc[cell_id]
            else:
                # 创建基本metadata
                cell_obs = pd.Series({'cell_id': cell_id})
            
            obs_list.append(cell_obs)
            
        except Exception as e:
            logger.warning(f"处理细胞 {cell_id} 失败: {e}")
            failed_cells.append(cell_id)
    
    if not X_list:
        raise ValueError("没有成功处理任何细胞")
    
    logger.info(f"成功处理 {len(X_list)} 个细胞，失败 {len(failed_cells)} 个")
    
    # 5. 构建矩阵
    X = np.vstack(X_list)
    obs_df_final = pd.DataFrame(obs_list)
    obs_df_final.index = [f"cell_{i}" for i in range(len(obs_df_final))]
    
    # 6. 创建 var DataFrame (特征名称)
    var_names = []
    for chr_label in chr_list:
        n_bins = chr_info[chr_label]
        n_tri = n_bins * (n_bins + 1) // 2
        for i in range(n_tri):
            var_names.append(f"{chr_label}_tri_{i}")
    
    var_df = pd.DataFrame(index=var_names)
    var_df['chromosome'] = [name.split('_')[0] for name in var_names]
    
    # 7. 创建 AnnData
    adata = ad.AnnData(X=X, obs=obs_df_final, var=var_df)
    
    # 添加一些基本信息
    adata.uns['stage'] = stage_dir.name
    adata.uns['n_chromosomes'] = len(chr_list)
    adata.uns['chromosomes'] = chr_list
    adata.uns['chr_bins'] = chr_info
    adata.uns['failed_cells'] = failed_cells
    
    # 8. 保存
    logger.info(f"保存到: {output_path}")
    logger.info(f"矩阵大小: {adata.shape}")
    logger.info(f"内存使用: {adata.X.nbytes / 1024**3:.2f} GB")
    
    adata.write_h5ad(output_path)
    logger.info("完成!")
    
    return adata


def main():
    parser = argparse.ArgumentParser(description='为单个 stage 构建 h5ad 文件')
    parser.add_argument('--stage-dir', required=True, help='stage 目录路径')
    parser.add_argument('--obs-xlsx', required=True, help='metadata Excel 文件')
    parser.add_argument('--output', required=True, help='输出 h5ad 文件路径')
    parser.add_argument('--max-cells', type=int, help='最大处理细胞数（用于测试）')
    
    args = parser.parse_args()
    
    stage_dir = Path(args.stage_dir)
    obs_xlsx = Path(args.obs_xlsx)
    output_path = Path(args.output)
    
    # 确保输出目录存在
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    build_stage_h5ad(stage_dir, obs_xlsx, output_path, args.max_cells)


if __name__ == '__main__':
    main()
