#!/usr/bin/env python3
"""
内存高效的 Hi-C NPZ 到 AnnData 转换器

使用 HDF5/AnnData 的 backed mode 实现单次遍历、分块写入，
避免将大量数据同时加载到内存中。

核心思想：
1. 预计算最终矩阵维度
2. 创建 file-backed AnnData 对象
3. 逐个细胞处理并直接写入硬盘
4. 内存峰值只等于单个细胞的数据量
"""

import numpy as np
import pandas as pd
import anndata as ad
import h5py
from pathlib import Path
import logging
from typing import List, Dict, Tuple, Optional
import warnings
import gc

# 导入现有的工具类
from hic_converter import (
    HiCNpzLoader, 
    UpperTriangleExtractor, 
    ChromosomeManager
)

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MemoryEfficientHiCConverter:
    """内存高效的 Hi-C 数据转换器"""
    
    def __init__(self, 
                 loader: HiCNpzLoader,
                 extractor: UpperTriangleExtractor,
                 chr_manager: ChromosomeManager,
                 file_pattern: str = "*_{}_pad1_std1.0_rp0.5_sqrtvc.npz"):
        """
        初始化转换器
        
        Parameters:
        -----------
        loader : HiCNpzLoader
            矩阵加载器
        extractor : UpperTriangleExtractor
            上三角提取器
        chr_manager : ChromosomeManager
            染色体管理器
        file_pattern : str
            文件名模式
        """
        self.loader = loader
        self.extractor = extractor
        self.chr_manager = chr_manager
        self.file_pattern = file_pattern
        
    def discover_all_cells(self, stage_dir: Path) -> List[Path]:
        """
        发现阶段目录中的所有细胞，按固定顺序排序
        
        Parameters:
        -----------
        stage_dir : Path
            阶段目录
            
        Returns:
        --------
        cell_dirs : list of Path
            按顺序排列的细胞目录列表
        """
        cell_dirs = []
        for item in stage_dir.iterdir():
            if item.is_dir():
                cell_dirs.append(item)
        
        # 确保固定顺序
        cell_dirs.sort(key=lambda x: x.name)
        logger.info(f"Found {len(cell_dirs)} cells in stage {stage_dir.name}")
        
        return cell_dirs
    
    def calculate_feature_dimensions(self, sample_cell_dir: Path) -> Tuple[int, List[str]]:
        """
        通过分析样本细胞计算特征维度
        
        Parameters:
        -----------
        sample_cell_dir : Path
            样本细胞目录
            
        Returns:
        --------
        n_features : int
            总特征数
        feature_names : list of str
            特征名称列表
        """
        logger.info(f"Calculating feature dimensions using sample cell: {sample_cell_dir.name}")
        
        # 验证染色体文件
        chr_files = self.chr_manager.validate_chromosome_files(sample_cell_dir, self.file_pattern)
        
        if not chr_files:
            raise ValueError(f"No chromosome files found in sample cell {sample_cell_dir}")
        
        total_features = 0
        feature_names = []
        
        # 按顺序处理每条染色体
        for chr_short in self.chr_manager.get_chromosome_order():
            chr_name = f"chr{chr_short}"
            if chr_name not in chr_files:
                logger.warning(f"Skipping missing chromosome {chr_name} in sample calculation")
                continue
            
            # 加载矩阵获取维度
            matrix = self.loader.load_sparse_matrix(chr_files[chr_name])
            n = matrix.shape[0]
            
            # 计算上三角元素数量
            if self.extractor.include_diagonal:
                chr_features = n * (n + 1) // 2
                indices = np.triu_indices(n, k=0)
            else:
                chr_features = n * (n - 1) // 2
                indices = np.triu_indices(n, k=1)
            
            total_features += chr_features
            
            # 生成特征名称
            chr_feature_names = [f"{chr_name}_bin{i}_bin{j}" 
                               for i, j in zip(indices[0], indices[1])]
            feature_names.extend(chr_feature_names)
            
            logger.info(f"Chromosome {chr_name}: {chr_features:,} features")
            
            # 释放内存
            del matrix
            gc.collect()
        
        logger.info(f"Total features calculated: {total_features:,}")
        logger.info(f"Estimated feature vector size: {total_features * 8 / 1024 / 1024:.2f} MB per cell")
        
        return total_features, feature_names
    
    def process_single_cell(self, cell_dir: Path) -> np.ndarray:
        """
        处理单个细胞，返回其特征向量
        
        Parameters:
        -----------
        cell_dir : Path
            细胞目录
            
        Returns:
        --------
        feature_vector : np.ndarray
            细胞的特征向量
        """
        # 验证染色体文件
        chr_files = self.chr_manager.validate_chromosome_files(cell_dir, self.file_pattern)
        
        if not chr_files:
            raise ValueError(f"No chromosome files found in {cell_dir}")
        
        cell_vectors = []
        
        # 按固定顺序处理每条染色体
        for chr_short in self.chr_manager.get_chromosome_order():
            chr_name = f"chr{chr_short}"
            if chr_name not in chr_files:
                logger.warning(f"Skipping missing chromosome {chr_name} for cell {cell_dir.name}")
                continue
            
            try:
                # 加载矩阵
                matrix = self.loader.load_sparse_matrix(chr_files[chr_name])
                
                # 使用稀疏方法提取上三角向量（内存效率更高）
                vector = self.extractor.extract_upper_triangle_sparse(matrix)
                cell_vectors.append(vector)
                
                # 立即释放矩阵内存
                del matrix
                
            except Exception as e:
                logger.error(f"Failed to process {chr_name} for cell {cell_dir.name}: {e}")
                raise
        
        if not cell_vectors:
            raise ValueError(f"No valid chromosome data found for cell {cell_dir.name}")
        
        # 水平拼接所有染色体的向量
        feature_vector = np.concatenate(cell_vectors)
        
        # 释放中间向量内存
        del cell_vectors
        gc.collect()
        
        return feature_vector
    
    def create_backed_anndata(self, 
                            output_file: Path,
                            n_cells: int,
                            n_features: int,
                            cell_ids: List[str],
                            stage_name: str) -> h5py.File:
        """
        直接创建 HDF5 文件用于存储数据（绕过 AnnData 的内存限制）
        
        Parameters:
        -----------
        output_file : Path
            输出文件路径
        n_cells : int
            细胞数量
        n_features : int
            特征数量
        cell_ids : list of str
            细胞ID列表
        stage_name : str
            发育阶段名称
            
        Returns:
        --------
        h5_file : h5py.File
            可写的 HDF5 文件对象
        """
        logger.info(f"Creating HDF5 file: {n_cells} cells × {n_features} features")
        
        # 创建 HDF5 文件
        h5_file = h5py.File(output_file, 'w')
        
        # 创建主数据矩阵（分块存储以提高性能）
        chunk_size = min(1000, n_features)
        h5_file.create_dataset(
            'X', 
            (n_cells, n_features), 
            dtype=np.float32, 
            chunks=(1, chunk_size),
            compression='gzip',
            compression_opts=1  # 轻度压缩以平衡速度和大小
        )
        
        # 创建观测信息
        obs_group = h5_file.create_group('obs')
        
        # 细胞ID
        cell_id_data = [id.encode('utf-8') for id in cell_ids]
        obs_group.create_dataset('_index', data=cell_id_data, dtype=h5py.string_dtype())
        obs_group.create_dataset('cell_id', data=cell_id_data, dtype=h5py.string_dtype())
        
        # 阶段信息
        stage_data = [stage_name.encode('utf-8')] * n_cells
        obs_group.create_dataset('stage', data=stage_data, dtype=h5py.string_dtype())
        
        # 添加元数据
        uns_group = h5_file.create_group('uns')
        processing_info = uns_group.create_group('processing_info')
        processing_info.attrs['conversion_method'] = 'memory_efficient_direct_hdf5'
        processing_info.attrs['include_diagonal'] = self.extractor.include_diagonal
        processing_info.attrs['stage'] = stage_name
        processing_info.attrs['n_cells'] = n_cells
        processing_info.attrs['n_features'] = n_features
        
        logger.info(f"HDF5 file structure created: {output_file}")
        
        return h5_file
    
    def convert_stage_memory_efficient(self, 
                                     stage_dir: Path,
                                     output_file: Path,
                                     max_cells: Optional[int] = None) -> Path:
        """
        内存高效地转换单个发育阶段的数据
        
        Parameters:
        -----------
        stage_dir : Path
            阶段目录
        output_file : Path
            输出文件路径
        max_cells : int, optional
            最大细胞数限制（用于测试）
            
        Returns:
        --------
        output_file : Path
            输出文件路径
        """
        stage_name = stage_dir.name
        logger.info(f"Starting memory-efficient conversion for stage: {stage_name}")
        
        # 第1步：发现所有细胞（固定顺序）
        all_cell_dirs = self.discover_all_cells(stage_dir)
        
        if max_cells is not None and len(all_cell_dirs) > max_cells:
            logger.info(f"Limiting to {max_cells} cells out of {len(all_cell_dirs)} for testing")
            all_cell_dirs = all_cell_dirs[:max_cells]
        
        if not all_cell_dirs:
            raise ValueError(f"No cells found in stage {stage_name}")
        
        n_cells = len(all_cell_dirs)
        cell_ids = [cell_dir.name for cell_dir in all_cell_dirs]
        
        # 第2步：计算特征维度（使用第一个细胞作为样本）
        n_features, _ = self.calculate_feature_dimensions(all_cell_dirs[0])
        
        # 第3步：创建 HDF5 文件
        logger.info("Creating HDF5 file...")
        h5_file = self.create_backed_anndata(
            output_file=output_file,
            n_cells=n_cells,
            n_features=n_features,
            cell_ids=cell_ids,
            stage_name=stage_name
        )
        
        # 第4步：逐个细胞处理并写入
        logger.info("Starting cell-by-cell processing...")
        
        try:
            for i, cell_dir in enumerate(all_cell_dirs):
                try:
                    logger.info(f"Processing cell {i+1}/{n_cells}: {cell_dir.name}")
                    
                    # 处理单个细胞
                    cell_features = self.process_single_cell(cell_dir)
                    
                    # 验证特征向量长度
                    if len(cell_features) != n_features:
                        logger.error(f"Feature dimension mismatch for cell {cell_dir.name}: "
                                   f"expected {n_features}, got {len(cell_features)}")
                        continue
                    
                    # 直接写入 HDF5 文件的相应行
                    h5_file['X'][i, :] = cell_features
                    
                    # 强制刷新到硬盘
                    h5_file.flush()
                    
                    # 释放内存
                    del cell_features
                    gc.collect()
                    
                    # 进度报告
                    if (i + 1) % 10 == 0:
                        memory_usage = self._get_memory_usage()
                        logger.info(f"Processed {i+1}/{n_cells} cells, Memory usage: {memory_usage:.2f} MB")
                    
                except Exception as e:
                    logger.error(f"Failed to process cell {cell_dir.name}: {e}")
                    continue
        
        finally:
            # 关闭 HDF5 文件
            logger.info("Finalizing and closing file...")
            h5_file.close()
        
        logger.info(f"Successfully converted stage {stage_name}")
        logger.info(f"Final output: {output_file}")
        logger.info(f"Final matrix shape: {n_cells} × {n_features}")
        
        return output_file
    
    def _get_memory_usage(self) -> float:
        """获取当前内存使用量（MB）"""
        import psutil
        import os
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
