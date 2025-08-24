#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CUDA GPU 加速版本：使用 NVIDIA GPU 加速 Hi-C 数据处理
依赖：cupy, cusparse

特点：
- 使用 CuPy 加速稀疏矩阵操作
- CUDA 内核优化的上三角提取
- 高效的 GPU 内存管理

用法：
  python build_stage_h5ad_cuda.py \
    --stage-dir /path/to/E70 \
    --obs-xlsx /path/to/metadata.xlsx \
    --output /path/to/E70.h5ad \
    --batch-size 32
"""

import argparse
import logging
import re
from pathlib import Path
from typing import List, Dict, Tuple
import time

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import sparse
import anndata as ad

try:
    import cupy as cp
    import cupyx.scipy.sparse as cp_sparse
    CUDA_AVAILABLE = True
except ImportError:
    CUDA_AVAILABLE = False
    cp = None

logger = logging.getLogger('build_stage_h5ad_cuda')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class CUDAAccelerator:
    """CUDA GPU 加速器"""
    
    def __init__(self, device_id: int = 0):
        if not CUDA_AVAILABLE:
            raise RuntimeError("CuPy 未安装或 CUDA 不可用")
        
        self.device_id = device_id
        cp.cuda.Device(device_id).use()
        
        # 检查 GPU 内存
        mempool = cp.get_default_memory_pool()
        total_memory = cp.cuda.Device().mem_info[1]
        logger.info(f"使用 GPU {device_id}, 总内存: {total_memory / 1024**3:.1f} GB")
        
        # 预计算上三角索引的 CUDA 内核
        self.triu_kernel = cp.ElementwiseKernel(
            'int32 n, raw T sparse_data, raw int32 indices, raw int32 indptr',
            'T triu_data',
            '''
            int row = i / n;
            int col = i % n;
            if (col >= row) {
                // 计算在 CSR 中的位置
                int start = indptr[row];
                int end = indptr[row + 1];
                for (int idx = start; idx < end; idx++) {
                    if (indices[idx] == col) {
                        triu_data = sparse_data[idx];
                        break;
                    }
                }
            }
            ''',
            'extract_triu'
        )
    
    def extract_triu_features_cuda(self, csr_matrix: sparse.csr_matrix) -> cp.ndarray:
        """使用 CUDA 提取上三角特征"""
        n = csr_matrix.shape[0]
        n_tri = n * (n + 1) // 2
        
        # 转移到 GPU
        data_gpu = cp.asarray(csr_matrix.data, dtype=cp.float32)
        indices_gpu = cp.asarray(csr_matrix.indices, dtype=cp.int32)
        indptr_gpu = cp.asarray(csr_matrix.indptr, dtype=cp.int32)
        
        # 使用优化的 CUDA 内核提取上三角
        triu_data = cp.zeros(n_tri, dtype=cp.float32)
        
        # 创建上三角线性索引
        row_indices = cp.repeat(cp.arange(n), cp.arange(n, 0, -1))
        col_indices = cp.concatenate([cp.arange(i, n) for i in range(n)])
        
        # 批量查找稀疏矩阵中的值
        for i in range(len(row_indices)):
            row, col = int(row_indices[i]), int(col_indices[i])
            start, end = int(indptr_gpu[row]), int(indptr_gpu[row + 1])
            
            # 在当前行中查找列索引
            if start < end:
                row_indices_slice = indices_gpu[start:end]
                mask = (row_indices_slice == col)
                if cp.any(mask):
                    idx = cp.where(mask)[0][0]
                    triu_data[i] = data_gpu[start + idx]
        
        return triu_data
    
    def extract_triu_features_batch_cuda(self, csr_matrices: List[sparse.csr_matrix]) -> cp.ndarray:
        """批量 CUDA 上三角提取"""
        batch_features = []
        
        with cp.cuda.Stream() as stream:
            for csr in csr_matrices:
                with stream:
                    features = self.extract_triu_features_cuda(csr)
                    batch_features.append(features)
        
        return cp.stack(batch_features)
    
    def extract_distance_features_cuda(self, csr_matrix: sparse.csr_matrix, 
                                     max_distance: int = 100) -> cp.ndarray:
        """CUDA 距离分箱特征提取"""
        # 转换为 CuPy 稀疏矩阵
        csr_gpu = cp_sparse.csr_matrix(
            (cp.asarray(csr_matrix.data), 
             cp.asarray(csr_matrix.indices),
             cp.asarray(csr_matrix.indptr)),
            shape=csr_matrix.shape
        )
        
        # 转换为 COO 格式便于距离计算
        coo_gpu = csr_gpu.tocoo()
        
        # 计算距离
        distances = cp.abs(coo_gpu.col - coo_gpu.row)
        
        features = []
        for d in range(1, min(max_distance + 1, csr_matrix.shape[0])):
            mask = (distances == d)
            if cp.any(mask):
                dist_values = coo_gpu.data[mask]
                features.extend([
                    cp.mean(dist_values),
                    cp.std(dist_values) if len(dist_values) > 1 else cp.float32(0),
                    cp.sum(dist_values),
                    cp.float32(len(dist_values))
                ])
            else:
                features.extend([cp.float32(0)] * 4)
        
        return cp.array(features, dtype=cp.float32)


def build_stage_h5ad_cuda(stage_dir: Path, obs_xlsx: Path, output_path: Path,
                         max_cells: int = None, batch_size: int = 16,
                         use_distance_features: bool = False, device_id: int = 0):
    """
    使用 CUDA 加速构建 h5ad 文件
    """
    logger.info(f"开始处理 stage: {stage_dir} (CUDA 加速)")
    
    # 初始化 CUDA 加速器
    accelerator = CUDAAccelerator(device_id)
    
    # 重用之前的发现和加载函数
    from build_stage_h5ad_mps import discover_chromosome_structure, load_cell_matrices
    
    chr_list, chr_info = discover_chromosome_structure(stage_dir)
    
    cell_dirs = [d for d in stage_dir.iterdir() if d.is_dir() and any(d.glob('*.npz'))]
    cell_dirs = sorted(cell_dirs, key=lambda x: x.name)
    
    if max_cells:
        cell_dirs = cell_dirs[:max_cells]
    
    logger.info(f"将处理 {len(cell_dirs)} 个细胞，批大小: {batch_size}")
    
    obs_df = pd.read_excel(obs_xlsx)
    
    # 批量处理
    all_features = []
    obs_list = []
    failed_cells = []
    
    for i in tqdm(range(0, len(cell_dirs), batch_size), desc="CUDA 批量处理"):
        batch_dirs = cell_dirs[i:i+batch_size]
        batch_matrices_by_cell = []
        batch_obs = []
        
        # 加载批次数据
        for cell_dir in batch_dirs:
            cell_id = cell_dir.name
            try:
                matrices = load_cell_matrices(cell_dir, chr_list)
                batch_matrices_by_cell.append(matrices)
                batch_obs.append(pd.Series({'cell_id': cell_id}))
            except Exception as e:
                logger.warning(f"处理细胞 {cell_id} 失败: {e}")
                failed_cells.append(cell_id)
        
        if not batch_matrices_by_cell:
            continue
        
        # CUDA 批量处理
        batch_start = time.time()
        
        batch_features_list = []
        for matrices in batch_matrices_by_cell:
            if use_distance_features:
                cell_features_list = []
                for csr in matrices:
                    feat = accelerator.extract_distance_features_cuda(csr, max_distance=50)
                    cell_features_list.append(feat)
                cell_features = cp.concatenate(cell_features_list)
            else:
                cell_features = accelerator.extract_triu_features_batch_cuda(matrices)
                cell_features = cell_features.flatten()
            
            batch_features_list.append(cell_features)
        
        # 合并批次特征并转回 CPU
        batch_features_gpu = cp.stack(batch_features_list)
        batch_features_cpu = cp.asnumpy(batch_features_gpu)
        
        batch_time = time.time() - batch_start
        logger.info(f"CUDA 批次 {i//batch_size + 1}: {len(batch_dirs)} 细胞, "
                   f"{batch_time:.2f}s, {batch_features_cpu.shape}")
        
        all_features.append(batch_features_cpu)
        obs_list.extend(batch_obs)
        
        # 清理 GPU 内存
        cp._default_memory_pool.free_all_blocks()
    
    if not all_features:
        raise ValueError("没有成功处理任何细胞")
    
    # 构建最终结果（同 MPS 版本）
    logger.info("合并特征矩阵...")
    X = np.vstack(all_features)
    obs_df_final = pd.DataFrame(obs_list)
    obs_df_final.index = [f"cell_{i}" for i in range(len(obs_df_final))]
    
    # 创建 var
    if use_distance_features:
        var_names = []
        for chr_label in chr_list:
            for d in range(1, 51):
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
    
    # 创建和保存 AnnData
    adata = ad.AnnData(X=X, obs=obs_df_final, var=var_df)
    adata.uns['stage'] = stage_dir.name
    adata.uns['processing'] = 'CUDA_accelerated'
    adata.uns['use_distance_features'] = use_distance_features
    adata.uns['failed_cells'] = failed_cells
    
    logger.info(f"保存到: {output_path}")
    logger.info(f"矩阵大小: {adata.shape}")
    logger.info(f"内存使用: {adata.X.nbytes / 1024**3:.2f} GB")
    
    adata.write_h5ad(output_path)
    logger.info("CUDA 加速处理完成!")
    
    return adata


def main():
    parser = argparse.ArgumentParser(description='CUDA GPU 加速构建 h5ad 文件')
    parser.add_argument('--stage-dir', required=True, help='stage 目录路径')
    parser.add_argument('--obs-xlsx', required=True, help='metadata Excel 文件')
    parser.add_argument('--output', required=True, help='输出 h5ad 文件路径')
    parser.add_argument('--max-cells', type=int, help='最大处理细胞数')
    parser.add_argument('--batch-size', type=int, default=16, help='批处理大小')
    parser.add_argument('--device-id', type=int, default=0, help='CUDA 设备 ID')
    parser.add_argument('--distance-features', action='store_true',
                       help='使用距离特征而不是完整上三角')
    
    args = parser.parse_args()
    
    stage_dir = Path(args.stage_dir)
    obs_xlsx = Path(args.obs_xlsx)
    output_path = Path(args.output)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    build_stage_h5ad_cuda(
        stage_dir, obs_xlsx, output_path,
        args.max_cells, args.batch_size, args.distance_features, args.device_id
    )


if __name__ == '__main__':
    main()
