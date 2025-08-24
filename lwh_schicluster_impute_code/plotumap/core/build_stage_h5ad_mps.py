#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MPS GPU 加速版本：使用 Apple Silicon GPU 加速 Hi-C 数据处理
依赖：pytorch (with MPS support)

特点：
- 使用 PyTorch MPS 后端加速稀疏矩阵操作
- 批量处理减少 CPU-GPU 数据传输
- 优化的内存管理

用法：
  python build_stage_h5ad_mps.py \
    --stage-dir /path/to/E70 \
    --obs-xlsx /path/to/metadata.xlsx \
    --output /path/to/E70.h5ad \
    --batch-size 16
"""

import argparse
import logging
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import time

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import sparse
import anndata as ad

# PyTorch for MPS acceleration
import torch
import torch.sparse

logger = logging.getLogger('build_stage_h5ad_mps')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class MPSAccelerator:
    """MPS GPU 加速器"""
    
    def __init__(self):
        if not torch.backends.mps.is_available():
            raise RuntimeError("MPS 不可用，请检查 PyTorch 安装和硬件支持")
        
        self.device = torch.device("mps")
        logger.info(f"使用设备: {self.device}")
        
        # 预计算上三角索引缓存
        self.triu_cache = {}
    
    def get_triu_indices(self, n: int) -> torch.Tensor:
        """获取或计算上三角索引"""
        if n not in self.triu_cache:
            # 在 GPU 上计算上三角索引
            rows, cols = torch.triu_indices(n, n, device=self.device)
            linear_idx = rows * n - rows * (rows - 1) // 2 + (cols - rows)
            self.triu_cache[n] = linear_idx
        return self.triu_cache[n]
    
    def extract_triu_features_batch(self, csr_matrices: List[sparse.csr_matrix]) -> torch.Tensor:
        """批量提取上三角特征"""
        batch_features = []
        
        for csr in csr_matrices:
            n = csr.shape[0]
            
            # 转换为 PyTorch 稀疏张量
            coo = csr.tocoo()
            indices = torch.from_numpy(np.vstack([coo.row, coo.col])).long()
            values = torch.from_numpy(coo.data).float()
            sparse_tensor = torch.sparse_coo_tensor(
                indices, values, csr.shape, device=self.device
            )
            
            # 提取上三角
            triu_indices = self.get_triu_indices(n)
            n_tri = n * (n + 1) // 2
            
            # 创建结果张量
            triu_features = torch.zeros(n_tri, device=self.device, dtype=torch.float32)
            
            # 获取上三角位置的值
            rows, cols = torch.triu_indices(n, n, device=self.device)
            if sparse_tensor._nnz() > 0:
                # 使用 sparse tensor 的 index_select 操作
                dense_triu = sparse_tensor.to_dense()[rows, cols]
                triu_features = dense_triu
            
            batch_features.append(triu_features)
        
        # 拼接所有特征
        return torch.cat(batch_features, dim=0)
    
    def extract_distance_features_batch(self, csr_matrices: List[sparse.csr_matrix], 
                                       max_distance: int = 100) -> torch.Tensor:
        """批量提取距离分箱特征（更快的替代方案）"""
        batch_features = []
        
        for csr in csr_matrices:
            features = []
            n = csr.shape[0]
            
            # 转换为 GPU 张量
            coo = csr.tocoo()
            indices = torch.from_numpy(np.vstack([coo.row, coo.col])).long().to(self.device)
            values = torch.from_numpy(coo.data).float().to(self.device)
            
            # 计算距离
            distances = torch.abs(indices[1] - indices[0])
            
            # 距离分箱特征
            for d in range(1, min(max_distance + 1, n)):
                mask = (distances == d)
                if mask.any():
                    dist_values = values[mask]
                    # 统计特征
                    features.extend([
                        dist_values.mean(),
                        dist_values.std() if len(dist_values) > 1 else torch.tensor(0.0, device=self.device),
                        dist_values.sum(),
                        torch.tensor(len(dist_values), device=self.device, dtype=torch.float32)
                    ])
                else:
                    features.extend([torch.tensor(0.0, device=self.device)] * 4)
            
            batch_features.append(torch.stack(features))
        
        return torch.stack(batch_features)


def extract_chr_label(filename: str) -> str:
    """从文件名提取染色体标签"""
    m = re.search(r"(chr(?:[0-9]+|X))", filename)
    return m.group(1) if m else None


def sort_chromosomes(chr_list: List[str]) -> List[str]:
    """对染色体进行排序"""
    def chr_key(label: str):
        if label == 'chrX':
            return 1000
        try:
            return int(label[3:])
        except:
            return 999
    return sorted(chr_list, key=chr_key)


def discover_chromosome_structure(stage_dir: Path) -> Tuple[List[str], Dict[str, int]]:
    """扫描染色体结构"""
    cell_dirs = [d for d in stage_dir.iterdir() if d.is_dir() and any(d.glob('*.npz'))]
    if not cell_dirs:
        raise ValueError(f"在 {stage_dir} 中未找到包含 .npz 文件的细胞目录")
    
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
    
    chr_list = sort_chromosomes(list(chr_info.keys()))
    
    logger.info(f"发现 {len(chr_list)} 条染色体: {chr_list}")
    return chr_list, chr_info


def load_cell_matrices(cell_dir: Path, chr_list: List[str]) -> List[sparse.csr_matrix]:
    """加载单个细胞的所有染色体矩阵"""
    matrices = []
    
    for chr_label in chr_list:
        npz_file = None
        for f in cell_dir.glob('*.npz'):
            if extract_chr_label(f.name) == chr_label:
                npz_file = f
                break
        
        if npz_file is None:
            raise FileNotFoundError(f"细胞 {cell_dir.name} 缺少染色体 {chr_label}")
        
        arr = np.load(npz_file, allow_pickle=True)
        csr = sparse.csr_matrix(
            (arr['data'], arr['indices'], arr['indptr']), 
            shape=tuple(arr['shape'])
        )
        matrices.append(csr)
    
    return matrices


def build_stage_h5ad_mps(stage_dir: Path, obs_xlsx: Path, output_path: Path, 
                        max_cells: int = None, batch_size: int = 8, 
                        use_distance_features: bool = False):
    """
    使用 MPS 加速构建 h5ad 文件
    """
    logger.info(f"开始处理 stage: {stage_dir} (MPS 加速)")
    
    # 初始化 MPS 加速器
    accelerator = MPSAccelerator()
    
    # 发现染色体结构
    chr_list, chr_info = discover_chromosome_structure(stage_dir)
    
    # 收集细胞目录
    cell_dirs = [d for d in stage_dir.iterdir() if d.is_dir() and any(d.glob('*.npz'))]
    cell_dirs = sorted(cell_dirs, key=lambda x: x.name)
    
    if max_cells:
        cell_dirs = cell_dirs[:max_cells]
    
    logger.info(f"将处理 {len(cell_dirs)} 个细胞，批大小: {batch_size}")
    
    # 加载metadata
    obs_df = pd.read_excel(obs_xlsx)
    
    # 批量处理细胞
    all_features = []
    obs_list = []
    failed_cells = []
    
    for i in tqdm(range(0, len(cell_dirs), batch_size), desc="批量处理"):
        batch_dirs = cell_dirs[i:i+batch_size]
        batch_matrices_by_cell = []
        batch_obs = []
        
        # 加载当前批次的所有矩阵
        for cell_dir in batch_dirs:
            cell_id = cell_dir.name
            try:
                matrices = load_cell_matrices(cell_dir, chr_list)
                batch_matrices_by_cell.append(matrices)
                
                # metadata
                cell_obs = pd.Series({'cell_id': cell_id})
                batch_obs.append(cell_obs)
                
            except Exception as e:
                logger.warning(f"处理细胞 {cell_id} 失败: {e}")
                failed_cells.append(cell_id)
        
        if not batch_matrices_by_cell:
            continue
        
        # GPU 批量处理
        batch_start = time.time()
        
        if use_distance_features:
            # 使用距离特征（更快）
            batch_features_list = []
            for matrices in batch_matrices_by_cell:
                cell_features = accelerator.extract_distance_features_batch(matrices, max_distance=50)
                cell_features_flat = cell_features.flatten()
                batch_features_list.append(cell_features_flat)
            
            batch_features = torch.stack(batch_features_list)
        else:
            # 使用完整上三角特征
            batch_features_list = []
            for matrices in batch_matrices_by_cell:
                cell_features = accelerator.extract_triu_features_batch(matrices)
                batch_features_list.append(cell_features)
            
            batch_features = torch.stack(batch_features_list)
        
        # 转回 CPU
        batch_features_cpu = batch_features.cpu().numpy()
        batch_time = time.time() - batch_start
        
        logger.info(f"批次 {i//batch_size + 1}: {len(batch_dirs)} 细胞, "
                   f"{batch_time:.2f}s, {batch_features_cpu.shape}")
        
        all_features.append(batch_features_cpu)
        obs_list.extend(batch_obs)
    
    if not all_features:
        raise ValueError("没有成功处理任何细胞")
    
    # 合并所有特征
    logger.info("合并特征矩阵...")
    X = np.vstack(all_features)
    obs_df_final = pd.DataFrame(obs_list)
    obs_df_final.index = [f"cell_{i}" for i in range(len(obs_df_final))]
    
    # 创建 var 信息
    if use_distance_features:
        var_names = []
        for chr_label in chr_list:
            for d in range(1, 51):  # max_distance=50
                for stat in ['mean', 'std', 'sum', 'count']:
                    var_names.append(f"{chr_label}_dist{d}_{stat}")
    else:
        var_names = []
        for chr_label in chr_list:
            n_bins = chr_info[chr_label]
            n_tri = n_bins * (n_bins + 1) // 2
            for i in range(n_tri):
                var_names.append(f"{chr_label}_tri_{i}")
    
    var_df = pd.DataFrame(index=var_names)
    var_df['chromosome'] = [name.split('_')[0] for name in var_names]
    
    # 创建 AnnData
    adata = ad.AnnData(X=X, obs=obs_df_final, var=var_df)
    adata.uns['stage'] = stage_dir.name
    adata.uns['processing'] = 'MPS_accelerated'
    adata.uns['use_distance_features'] = use_distance_features
    adata.uns['failed_cells'] = failed_cells
    
    # 保存
    logger.info(f"保存到: {output_path}")
    logger.info(f"矩阵大小: {adata.shape}")
    logger.info(f"内存使用: {adata.X.nbytes / 1024**3:.2f} GB")
    
    adata.write_h5ad(output_path)
    logger.info("完成!")
    
    return adata


def main():
    parser = argparse.ArgumentParser(description='MPS GPU 加速构建 h5ad 文件')
    parser.add_argument('--stage-dir', required=True, help='stage 目录路径')
    parser.add_argument('--obs-xlsx', required=True, help='metadata Excel 文件')
    parser.add_argument('--output', required=True, help='输出 h5ad 文件路径')
    parser.add_argument('--max-cells', type=int, help='最大处理细胞数')
    parser.add_argument('--batch-size', type=int, default=8, help='批处理大小')
    parser.add_argument('--distance-features', action='store_true', 
                       help='使用距离特征而不是完整上三角（更快但信息量略少）')
    
    args = parser.parse_args()
    
    stage_dir = Path(args.stage_dir)
    obs_xlsx = Path(args.obs_xlsx)
    output_path = Path(args.output)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    build_stage_h5ad_mps(
        stage_dir, obs_xlsx, output_path, 
        args.max_cells, args.batch_size, args.distance_features
    )


if __name__ == '__main__':
    main()
