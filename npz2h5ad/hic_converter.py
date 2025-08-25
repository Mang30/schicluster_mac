#!/usr/bin/env python3
"""
Hi-C NPZ 数据转换为 AnnData 格式的核心工具库

此模块提供了将单细胞 Hi-C 数据从 NPZ 格式转换为 AnnData 格式的核心功能。
主要功能包括：
1. 加载稀疏矩阵数据
2. 提取上三角矩阵
3. 生成特征向量
4. 构建 AnnData 对象
"""

import numpy as np
import scipy.sparse as sp
import pandas as pd
import anndata as ad
from pathlib import Path
import logging
from typing import List, Dict, Tuple, Optional, Union
import warnings

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class HiCNpzLoader:
    """Hi-C NPZ 文件加载器"""
    
    def __init__(self):
        pass
    
    def load_sparse_matrix(self, filepath: Union[str, Path]) -> sp.spmatrix:
        """
        从 NPZ 文件加载稀疏矩阵
        
        Parameters:
        -----------
        filepath : str or Path
            NPZ 文件路径
            
        Returns:
        --------
        matrix : scipy.sparse matrix
            加载的稀疏矩阵
        """
        try:
            data = np.load(filepath, allow_pickle=True)
            
            format_type = data['format'].item().decode('utf-8')
            shape = tuple(data['shape'])
            matrix_data = data['data']
            
            if format_type == 'coo':
                row = data['row']
                col = data['col']
                matrix = sp.coo_matrix((matrix_data, (row, col)), shape=shape)
            elif format_type == 'csr':
                indices = data['indices']
                indptr = data['indptr']
                matrix = sp.csr_matrix((matrix_data, indices, indptr), shape=shape)
            else:
                raise ValueError(f"Unsupported sparse matrix format: {format_type}")
            
            data.close()
            return matrix
            
        except Exception as e:
            logger.error(f"Failed to load matrix from {filepath}: {e}")
            raise

class UpperTriangleExtractor:
    """上三角矩阵提取器"""
    
    def __init__(self, include_diagonal: bool = True):
        """
        初始化上三角提取器
        
        Parameters:
        -----------
        include_diagonal : bool, default=True
            是否包含对角线元素
        """
        self.include_diagonal = include_diagonal
    
    def extract_upper_triangle(self, matrix: sp.spmatrix) -> np.ndarray:
        """
        提取矩阵的上三角部分并展平为一维向量
        
        Parameters:
        -----------
        matrix : scipy.sparse matrix
            输入的方形矩阵
            
        Returns:
        --------
        vector : np.ndarray
            展平的上三角向量
        """
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError(f"Matrix must be square, got shape {matrix.shape}")
        
        n = matrix.shape[0]
        
        # 转换为密集矩阵（对于大矩阵可能需要优化）
        if sp.issparse(matrix):
            dense_matrix = matrix.toarray()
        else:
            dense_matrix = matrix
        
        # 获取上三角索引
        if self.include_diagonal:
            row_indices, col_indices = np.triu_indices(n, k=0)
        else:
            row_indices, col_indices = np.triu_indices(n, k=1)
        
        # 提取上三角元素
        upper_tri_vector = dense_matrix[row_indices, col_indices]
        
        return upper_tri_vector
    
    def extract_upper_triangle_sparse(self, matrix: sp.spmatrix) -> np.ndarray:
        """
        高效提取稀疏矩阵的上三角部分（适用于大矩阵）
        
        Parameters:
        -----------
        matrix : scipy.sparse matrix
            输入的稀疏方形矩阵
            
        Returns:
        --------
        vector : np.ndarray
            展平的上三角向量
        """
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError(f"Matrix must be square, got shape {matrix.shape}")
        
        n = matrix.shape[0]
        
        # 转换为COO格式以便操作
        if not sp.isspmatrix_coo(matrix):
            matrix_coo = matrix.tocoo()
        else:
            matrix_coo = matrix
        
        # 筛选上三角元素
        if self.include_diagonal:
            mask = matrix_coo.row <= matrix_coo.col
        else:
            mask = matrix_coo.row < matrix_coo.col
        
        # 创建上三角稀疏矩阵
        upper_tri_matrix = sp.coo_matrix(
            (matrix_coo.data[mask], 
             (matrix_coo.row[mask], matrix_coo.col[mask])), 
            shape=matrix.shape
        )
        
        # 转换为密集向量
        if self.include_diagonal:
            row_indices, col_indices = np.triu_indices(n, k=0)
        else:
            row_indices, col_indices = np.triu_indices(n, k=1)
        
        # 创建结果向量
        result_vector = np.zeros(len(row_indices), dtype=matrix.dtype)
        
        # 填充非零值
        upper_tri_dense = upper_tri_matrix.toarray()
        result_vector = upper_tri_dense[row_indices, col_indices]
        
        return result_vector

class ChromosomeManager:
    """染色体管理器"""
    
    def __init__(self, chromosomes: Optional[List[str]] = None):
        """
        初始化染色体管理器
        
        Parameters:
        -----------
        chromosomes : list of str, optional
            染色体列表，如果为None则使用默认顺序
        """
        if chromosomes is None:
            # 默认染色体顺序（不包含 "chr" 前缀，在匹配时添加）
            self.chromosomes = [
                '1', '2', '3', '4', '5', '6', '7', '8',
                '9', '10', '11', '12', '13', '14', '15',
                '16', '17', '18', '19', 'X'
            ]
        else:
            self.chromosomes = chromosomes
    
    def get_chromosome_order(self) -> List[str]:
        """获取染色体顺序"""
        return self.chromosomes.copy()
    
    def validate_chromosome_files(self, cell_dir: Path, 
                                file_pattern: str = "*_{}_pad1_std1.0_rp0.5_sqrtvc.npz") -> Dict[str, Path]:
        """
        验证细胞目录中的染色体文件
        
        Parameters:
        -----------
        cell_dir : Path
            细胞目录
        file_pattern : str
            文件名模式，{} 会被染色体名替换
            
        Returns:
        --------
        chr_files : dict
            染色体名到文件路径的映射
        """
        chr_files = {}
        missing_chrs = []
        
        for chr_name in self.chromosomes:
            # 构建完整的染色体名称（添加 chr 前缀）
            full_chr_name = f"chr{chr_name}"
            pattern = file_pattern.format(full_chr_name)
            files = list(cell_dir.glob(pattern))
            
            if len(files) == 0:
                missing_chrs.append(full_chr_name)
            elif len(files) == 1:
                chr_files[full_chr_name] = files[0]
            else:
                # 如果有多个匹配文件，选择第一个并发出警告
                logger.warning(f"Multiple files found for {full_chr_name} in {cell_dir}, using {files[0]}")
                chr_files[full_chr_name] = files[0]
        
        if missing_chrs:
            logger.warning(f"Missing chromosome files in {cell_dir}: {missing_chrs}")
        
        return chr_files

class FeatureGenerator:
    """特征向量生成器"""
    
    def __init__(self, 
                 loader: HiCNpzLoader,
                 extractor: UpperTriangleExtractor,
                 chr_manager: ChromosomeManager):
        """
        初始化特征生成器
        
        Parameters:
        -----------
        loader : HiCNpzLoader
            矩阵加载器
        extractor : UpperTriangleExtractor
            上三角提取器
        chr_manager : ChromosomeManager
            染色体管理器
        """
        self.loader = loader
        self.extractor = extractor
        self.chr_manager = chr_manager
    
    def generate_cell_features(self, cell_dir: Path, 
                             file_pattern: str = "*_{}_pad1_std1.0_rp0.5_sqrtvc.npz") -> Tuple[np.ndarray, List[str]]:
        """
        为单个细胞生成特征向量
        
        Parameters:
        -----------
        cell_dir : Path
            细胞目录
        file_pattern : str
            文件名模式
            
        Returns:
        --------
        feature_vector : np.ndarray
            细胞的特征向量
        feature_info : list of str
            特征信息列表
        """
        # 验证染色体文件
        chr_files = self.chr_manager.validate_chromosome_files(cell_dir, file_pattern)
        
        if not chr_files:
            raise ValueError(f"No chromosome files found in {cell_dir}")
        
        cell_vectors = []
        feature_info = []
        
        # 按顺序处理每条染色体
        for chr_short in self.chr_manager.get_chromosome_order():
            chr_name = f"chr{chr_short}"  # 构建完整染色体名
            if chr_name not in chr_files:
                logger.warning(f"Skipping missing chromosome {chr_name} for cell {cell_dir.name}")
                continue
            
            try:
                # 加载矩阵
                matrix = self.loader.load_sparse_matrix(chr_files[chr_name])
                
                # 提取上三角向量
                if matrix.nnz > 100000:  # 对于大矩阵使用稀疏方法
                    vector = self.extractor.extract_upper_triangle_sparse(matrix)
                else:
                    vector = self.extractor.extract_upper_triangle(matrix)
                
                cell_vectors.append(vector)
                
                # 生成特征信息
                n = matrix.shape[0]
                if self.extractor.include_diagonal:
                    indices = np.triu_indices(n, k=0)
                else:
                    indices = np.triu_indices(n, k=1)
                
                chr_feature_info = [f"{chr_name}_bin{i}_bin{j}" 
                                  for i, j in zip(indices[0], indices[1])]
                feature_info.extend(chr_feature_info)
                
                logger.info(f"Processed {chr_name}: shape {matrix.shape}, "
                          f"features {len(vector)}")
                
            except Exception as e:
                logger.error(f"Failed to process {chr_name} for cell {cell_dir.name}: {e}")
                raise
        
        if not cell_vectors:
            raise ValueError(f"No valid chromosome data found for cell {cell_dir.name}")
        
        # 拼接所有染色体的向量
        feature_vector = np.concatenate(cell_vectors)
        
        return feature_vector, feature_info

class AnnDataBuilder:
    """AnnData 对象构建器"""
    
    def __init__(self):
        pass
    
    def build_feature_dataframe(self, feature_info: List[str]) -> pd.DataFrame:
        """
        构建特征信息 DataFrame
        
        Parameters:
        -----------
        feature_info : list of str
            特征信息列表，格式为 "chr_binI_binJ"
            
        Returns:
        --------
        var_df : pd.DataFrame
            特征信息 DataFrame
        """
        # 解析特征信息
        chromosomes = []
        bin1_indices = []
        bin2_indices = []
        
        for feature in feature_info:
            parts = feature.split('_')
            if len(parts) >= 3:
                chr_name = parts[0]
                bin1 = int(parts[1].replace('bin', ''))
                bin2 = int(parts[2].replace('bin', ''))
                
                chromosomes.append(chr_name)
                bin1_indices.append(bin1)
                bin2_indices.append(bin2)
            else:
                logger.warning(f"Invalid feature format: {feature}")
                chromosomes.append('unknown')
                bin1_indices.append(-1)
                bin2_indices.append(-1)
        
        var_df = pd.DataFrame({
            'chromosome': chromosomes,
            'bin1': bin1_indices,
            'bin2': bin2_indices,
            'feature_name': feature_info
        })
        
        var_df.index = feature_info
        
        return var_df
    
    def build_observation_dataframe(self, cell_ids: List[str], 
                                  stage_info: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """
        构建观测信息 DataFrame
        
        Parameters:
        -----------
        cell_ids : list of str
            细胞ID列表
        stage_info : dict, optional
            细胞阶段信息
            
        Returns:
        --------
        obs_df : pd.DataFrame
            观测信息 DataFrame
        """
        obs_data = {'cell_id': cell_ids}
        
        if stage_info:
            obs_data['stage'] = [stage_info.get(cell_id, 'unknown') for cell_id in cell_ids]
        
        obs_df = pd.DataFrame(obs_data)
        obs_df.index = cell_ids
        
        return obs_df
    
    def create_anndata(self, 
                      feature_matrix: np.ndarray,
                      cell_ids: List[str],
                      feature_info: List[str],
                      stage_info: Optional[Dict[str, str]] = None) -> ad.AnnData:
        """
        创建 AnnData 对象
        
        Parameters:
        -----------
        feature_matrix : np.ndarray
            特征矩阵 (n_cells x n_features)
        cell_ids : list of str
            细胞ID列表
        feature_info : list of str
            特征信息列表
        stage_info : dict, optional
            细胞阶段信息
            
        Returns:
        --------
        adata : anndata.AnnData
            AnnData 对象
        """
        # 构建观测和特征信息
        obs_df = self.build_observation_dataframe(cell_ids, stage_info)
        var_df = self.build_feature_dataframe(feature_info)
        
        # 创建 AnnData 对象
        adata = ad.AnnData(
            X=feature_matrix,
            obs=obs_df,
            var=var_df
        )
        
        # 添加元数据
        adata.uns['processing_info'] = {
            'conversion_method': 'upper_triangle_concatenation',
            'include_diagonal': True,  # 这个应该从 extractor 获取
            'n_chromosomes': len(set(var_df['chromosome'])),
            'chromosomes': sorted(list(set(var_df['chromosome'])))
        }
        
        logger.info(f"Created AnnData object: {adata.shape}")
        logger.info(f"Chromosomes: {adata.uns['processing_info']['chromosomes']}")
        
        return adata
