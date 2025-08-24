#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hi-C特化的AnnData构建模块
从合并后的Hi-C矩阵构建适合单细胞分析的AnnData对象
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import load_npz, csr_matrix
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import json
import time

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class HiCAnnDataBuilder:
    """Hi-C特化的AnnData构建器"""
    
    def __init__(self, chrom_sizes_file: str, resolution: int = 40000):
        """
        初始化构建器
        
        Parameters:
        -----------
        chrom_sizes_file : str
            染色体大小文件路径
        resolution : int
            Hi-C数据分辨率，默认100kb
        """
        self.resolution = resolution
        self.chrom_sizes = self._load_chrom_sizes(chrom_sizes_file)
        self.chrom_order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
                           'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                           'chr16', 'chr17', 'chr18', 'chr19', 'chrX']
        
        # 计算bin信息
        self.bin_info = self._create_bin_info()
        self.total_bins = len(self.bin_info)
        
        logger.info(f"初始化完成: 分辨率 {resolution}, 总bins {self.total_bins}")
        
    def _load_chrom_sizes(self, chrom_sizes_file: str) -> Dict[str, int]:
        """加载染色体大小信息"""
        chrom_sizes = {}
        try:
            df = pd.read_csv(chrom_sizes_file, sep='\t', header=None, names=['chrom', 'size'])
            for _, row in df.iterrows():
                chrom_sizes[row['chrom']] = int(row['size'])
            logger.info(f"加载染色体大小信息: {len(chrom_sizes)} 个染色体")
            return chrom_sizes
        except Exception as e:
            logger.error(f"加载染色体大小文件失败: {e}")
            raise
            
    def _create_bin_info(self) -> pd.DataFrame:
        """创建bin信息表"""
        bin_data = []
        bin_id = 0
        
        for chrom in self.chrom_order:
            if chrom not in self.chrom_sizes:
                continue
                
            chrom_size = self.chrom_sizes[chrom]
            n_bins = chrom_size // self.resolution + 1
            
            for i in range(n_bins):
                start_pos = i * self.resolution
                end_pos = min((i + 1) * self.resolution, chrom_size)
                
                bin_data.append({
                    'bin_id': bin_id,
                    'chrom': chrom,
                    'start': start_pos,
                    'end': end_pos,
                    'size': end_pos - start_pos
                })
                bin_id += 1
                
        bin_df = pd.DataFrame(bin_data)
        logger.info(f"创建bin信息: {len(bin_df)} 个bins")
        return bin_df
        
    def load_cell_metadata(self, stage_files_mapping: str, 
                          metadata_file: Optional[str] = None) -> pd.DataFrame:
        """
        加载细胞metadata信息
        
        Parameters:
        -----------
        stage_files_mapping : str
            stage文件映射CSV路径
        metadata_file : str, optional
            额外的metadata文件路径
            
        Returns:
        --------
        pd.DataFrame : 细胞metadata信息
        """
        try:
            # 加载stage mapping文件
            mapping_df = pd.read_csv(stage_files_mapping)
            
            # 从文件名中提取细胞名称
            mapping_df['cell_name'] = mapping_df['cellname']
            
            # 创建基本的metadata
            metadata_cols = ['cell_name', 'stage', 'celltype']
            if 'celltype' not in mapping_df.columns:
                mapping_df['celltype'] = 'Unknown'
                
            metadata_df = mapping_df[metadata_cols].drop_duplicates()
            
            # 如果提供了额外的metadata文件，进行合并
            if metadata_file and os.path.exists(metadata_file):
                try:
                    if metadata_file.endswith('.xlsx'):
                        extra_metadata = pd.read_excel(metadata_file)
                    else:
                        extra_metadata = pd.read_csv(metadata_file)
                    
                    # 尝试按细胞名称合并
                    merge_key = None
                    for col in ['Cellname', 'cell_name', 'cellname']:
                        if col in extra_metadata.columns:
                            merge_key = col
                            break
                            
                    if merge_key:
                        metadata_df = pd.merge(metadata_df, extra_metadata, 
                                             left_on='cell_name', right_on=merge_key, how='left')
                        logger.info(f"合并了额外的metadata信息: {extra_metadata.shape}")
                except Exception as e:
                    logger.warning(f"合并额外metadata失败，继续使用基本信息: {e}")
            
            logger.info(f"加载细胞metadata: {len(metadata_df)} 个细胞")
            return metadata_df
            
        except Exception as e:
            logger.error(f"加载细胞metadata失败: {e}")
            raise
            
    def build_stage_anndata(self, merged_files: Dict[str, str], stage_name: str,
                           metadata_df: pd.DataFrame, output_path: str,
                           normalization_method: str = 'log_cpm') -> bool:
        """
        为单个stage构建AnnData对象
        
        Parameters:
        -----------
        merged_files : Dict[str, str]
            细胞名称到合并矩阵文件路径的映射
        stage_name : str
            stage名称
        metadata_df : pd.DataFrame
            细胞metadata
        output_path : str
            输出h5ad文件路径
        normalization_method : str
            标准化方法：'log_cpm', 'zscore', 'quantile'
            
        Returns:
        --------
        bool : 是否成功构建
        """
        logger.info(f"开始构建Stage {stage_name} 的AnnData对象")
        start_time = time.time()
        
        try:
            # 检查输出文件是否已存在
            if os.path.exists(output_path):
                logger.info(f"AnnData文件已存在，跳过: {output_path}")
                return True
                
            # 加载所有细胞的矩阵
            cell_names = list(merged_files.keys())
            matrices = []
            valid_cells = []
            
            logger.info(f"加载 {len(cell_names)} 个细胞的矩阵数据")
            
            for cell_name in cell_names:
                try:
                    matrix_path = merged_files[cell_name]
                    matrix = load_npz(matrix_path)
                    
                    # 检查矩阵尺寸
                    if matrix.shape[0] != self.total_bins or matrix.shape[1] != self.total_bins:
                        logger.warning(f"细胞 {cell_name} 矩阵尺寸不匹配: {matrix.shape}")
                        continue
                        
                    matrices.append(matrix)
                    valid_cells.append(cell_name)
                    
                    if len(valid_cells) % 10 == 0:
                        logger.info(f"已加载 {len(valid_cells)} 个细胞矩阵")
                        
                except Exception as e:
                    logger.warning(f"加载细胞 {cell_name} 矩阵失败: {e}")
                    continue
                    
            if not matrices:
                logger.error(f"Stage {stage_name} 没有有效的细胞矩阵")
                return False
                
            logger.info(f"成功加载 {len(matrices)} 个细胞矩阵")
            
            # 将Hi-C矩阵转换为特征向量（上三角部分）
            features_list = []
            
            logger.info("提取Hi-C特征向量")
            for i, matrix in enumerate(matrices):
                # 提取上三角部分作为特征
                features = self._extract_hic_features(matrix)
                features_list.append(features)
                
                if (i + 1) % 10 == 0:
                    logger.info(f"已处理 {i + 1}/{len(matrices)} 个矩阵")
                    
            # 合并特征矩阵
            X = np.vstack(features_list)
            logger.info(f"特征矩阵形状: {X.shape}")
            
            # 创建obs (细胞信息)
            obs_df = metadata_df[metadata_df['cell_name'].isin(valid_cells)].copy()
            obs_df = obs_df.set_index('cell_name')
            obs_df = obs_df.loc[valid_cells]  # 保证顺序一致
            
            # 创建var (特征信息)
            var_df = self._create_var_info()
            
            # 创建AnnData对象
            adata = sc.AnnData(X=X, obs=obs_df, var=var_df)
            
            # 添加基本信息
            adata.uns['stage'] = stage_name
            adata.uns['resolution'] = self.resolution
            adata.uns['n_chromosomes'] = len(self.chrom_order)
            adata.uns['total_bins'] = self.total_bins
            
            # 数据标准化
            self._normalize_data(adata, method=normalization_method)
            
            # 保存AnnData对象
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            adata.write(output_path)
            
            elapsed_time = time.time() - start_time
            logger.info(f"Stage {stage_name} AnnData构建完成，耗时: {elapsed_time:.2f}秒")
            logger.info(f"AnnData信息: {adata.n_obs} 细胞, {adata.n_vars} 特征")
            
            return True
            
        except Exception as e:
            logger.error(f"构建Stage {stage_name} AnnData失败: {e}")
            return False
            
    def _extract_hic_features(self, matrix: csr_matrix) -> np.ndarray:
        """
        从Hi-C矩阵提取特征向量
        
        Parameters:
        -----------
        matrix : csr_matrix
            Hi-C接触矩阵
            
        Returns:
        --------
        np.ndarray : 特征向量
        """
        # 方法1: 提取上三角部分
        n = matrix.shape[0]
        triu_indices = np.triu_indices(n, k=1)  # k=1排除对角线
        features = matrix[triu_indices].toarray().flatten()
        
        # 可选：添加对角线特征
        diag_features = matrix.diagonal()
        
        # 可选：添加统计特征
        row_sums = np.array(matrix.sum(axis=1)).flatten()
        col_sums = np.array(matrix.sum(axis=0)).flatten()
        
        # 合并特征（这里只使用上三角部分）
        all_features = features
        
        return all_features
        
    def _create_var_info(self) -> pd.DataFrame:
        """创建特征信息表"""
        # 为上三角部分的每个元素创建特征信息
        n = self.total_bins
        var_data = []
        
        # 上三角索引
        triu_indices = np.triu_indices(n, k=1)
        
        for i, (row_idx, col_idx) in enumerate(zip(triu_indices[0], triu_indices[1])):
            # 获取对应的bin信息
            bin1_info = self.bin_info.iloc[row_idx]
            bin2_info = self.bin_info.iloc[col_idx]
            
            # 创建特征名称
            feature_name = f"{bin1_info['chrom']}:{bin1_info['start']}-{bin1_info['end']}_" \
                          f"{bin2_info['chrom']}:{bin2_info['start']}-{bin2_info['end']}"
            
            # 计算距离（如果是同一染色体）
            if bin1_info['chrom'] == bin2_info['chrom']:
                distance = abs(bin2_info['start'] - bin1_info['start'])
                interaction_type = 'intra'
            else:
                distance = np.inf
                interaction_type = 'inter'
                
            var_data.append({
                'feature_name': feature_name,
                'bin1_chrom': bin1_info['chrom'],
                'bin1_start': bin1_info['start'],
                'bin1_end': bin1_info['end'],
                'bin2_chrom': bin2_info['chrom'],
                'bin2_start': bin2_info['start'],
                'bin2_end': bin2_info['end'],
                'distance': distance,
                'interaction_type': interaction_type
            })
            
        var_df = pd.DataFrame(var_data)
        var_df.index = var_df['feature_name']
        
        logger.info(f"创建特征信息: {len(var_df)} 个特征")
        return var_df
        
    def _normalize_data(self, adata: sc.AnnData, method: str = 'log_cpm'):
        """
        标准化Hi-C数据
        
        Parameters:
        -----------
        adata : sc.AnnData
            AnnData对象
        method : str
            标准化方法
        """
        logger.info(f"使用方法 {method} 标准化数据")
        
        if method == 'log_cpm':
            # CPM标准化 + log1p转换
            sc.pp.normalize_total(adata, target_sum=1e6)  # CPM
            sc.pp.log1p(adata)
            
        elif method == 'zscore':
            # Z-score标准化
            from sklearn.preprocessing import StandardScaler
            scaler = StandardScaler()
            adata.X = scaler.fit_transform(adata.X)
            
        elif method == 'quantile':
            # 分位数标准化
            from sklearn.preprocessing import QuantileTransformer
            transformer = QuantileTransformer(n_quantiles=1000, random_state=42)
            adata.X = transformer.fit_transform(adata.X)
            
        else:
            logger.warning(f"未知的标准化方法: {method}，跳过标准化")
            
        logger.info("数据标准化完成")
        

def main():
    """测试主函数"""
    # 配置路径
    chrom_sizes_file = "/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes.txt"
    merged_base_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/merged_matrices"
    stage_files_mapping = "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/plotumap/config/stage_files_mapping.csv"
    output_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/hic_h5ad_files"
    
    # 创建构建器
    builder = HiCAnnDataBuilder(chrom_sizes_file, resolution=40000)
    
    # 加载metadata
    metadata_df = builder.load_cell_metadata(stage_files_mapping)
    
    # 测试阶段
    test_stages = ['E70', 'E80', 'EX05']
    
    for stage in test_stages:
        stage_dir = os.path.join(merged_base_dir, stage)
        if os.path.exists(stage_dir):
            # 收集合并后的矩阵文件
            merged_files = {}
            for f in os.listdir(stage_dir):
                if f.endswith('_merged.npz'):
                    cell_name = f.replace('_merged.npz', '')
                    merged_files[cell_name] = os.path.join(stage_dir, f)
                    
            if merged_files:
                output_path = os.path.join(output_dir, f"{stage}_hic.h5ad")
                success = builder.build_stage_anndata(merged_files, stage, metadata_df, output_path)
                logger.info(f"Stage {stage} AnnData构建{'成功' if success else '失败'}")
            else:
                logger.warning(f"Stage {stage} 没有找到合并后的矩阵文件")
        else:
            logger.warning(f"Stage目录不存在: {stage_dir}")


if __name__ == "__main__":
    main()