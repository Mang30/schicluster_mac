#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
鲁棒的Hi-C矩阵合并模块
改进内存管理和错误处理
"""

import os
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, save_npz, load_npz, block_diag
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import gc
import time

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RobustHiCMatrixMerger:
    """鲁棒的Hi-C矩阵合并器，改进内存管理和错误处理"""
    
    def __init__(self, chrom_sizes_file: str, resolution: int = 40000):
        """
        初始化合并器
        
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
        
        # 计算每个染色体的bin数和偏移量
        self.chrom_bins = {}
        self.chrom_offsets = {}
        current_offset = 0
        
        for chrom in self.chrom_order:
            if chrom in self.chrom_sizes:
                bins = self.chrom_sizes[chrom] // resolution + 1
                self.chrom_bins[chrom] = bins
                self.chrom_offsets[chrom] = current_offset
                current_offset += bins
            
        self.total_bins = current_offset
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
            
    def _load_single_chromosome(self, file_path: str, chrom: str) -> Optional[csr_matrix]:
        """安全加载单个染色体矩阵"""
        try:
            if not os.path.exists(file_path):
                logger.warning(f"染色体文件不存在: {file_path}")
                return None
                
            # 检查文件大小
            file_size = os.path.getsize(file_path)
            if file_size == 0:
                logger.warning(f"染色体文件为空: {file_path}")
                return None
                
            # 加载矩阵
            matrix = load_npz(file_path)
            
            # 验证矩阵尺寸
            expected_bins = self.chrom_bins.get(chrom, 0)
            if matrix.shape[0] != expected_bins or matrix.shape[1] != expected_bins:
                logger.warning(f"染色体 {chrom} 矩阵尺寸不匹配: {matrix.shape} vs 期望 ({expected_bins}, {expected_bins})")
                # 尝试调整尺寸
                if matrix.shape[0] > expected_bins and matrix.shape[1] > expected_bins:
                    matrix = matrix[:expected_bins, :expected_bins]
                    logger.info(f"  已截断矩阵到期望尺寸")
                else:
                    # 如果太小，填充零矩阵
                    from scipy.sparse import coo_matrix
                    new_matrix = coo_matrix((expected_bins, expected_bins))
                    min_size = min(matrix.shape[0], expected_bins)
                    new_matrix = new_matrix.tocsr()
                    new_matrix[:min_size, :min_size] = matrix[:min_size, :min_size]
                    matrix = new_matrix
                    logger.info(f"  已填充矩阵到期望尺寸")
                    
            logger.debug(f"成功加载 {chrom}: {matrix.shape}, 非零元素: {matrix.nnz:,}")
            return matrix
            
        except Exception as e:
            logger.error(f"加载染色体 {chrom} 失败: {e}")
            return None
            
    def merge_cell_chromosomes_robust(self, cell_dir: str, cell_name: str, 
                                    output_path: str, use_memory_efficient: bool = True) -> bool:
        """
        鲁棒的细胞染色体矩阵合并
        
        Parameters:
        -----------
        cell_dir : str
            细胞目录路径
        cell_name : str
            细胞名称
        output_path : str
            输出合并矩阵的路径
        use_memory_efficient : bool
            是否使用内存高效模式
            
        Returns:
        --------
        bool : 是否成功合并
        """
        logger.info(f"开始合并细胞 {cell_name} 的染色体矩阵（鲁棒模式）")
        start_time = time.time()
        
        try:
            # 检查输出文件是否已存在
            if os.path.exists(output_path):
                logger.info(f"合并矩阵已存在，跳过: {output_path}")
                return True
                
            # 创建输出目录
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # 收集可用的染色体文件
            available_chroms = []
            chrom_matrices = []
            
            for chrom in self.chrom_order:
                if chrom not in self.chrom_sizes:
                    continue
                    
                # 构造文件路径
                file_pattern = f"{cell_name}_{chrom}_pad1_std1.0_rp0.5_sqrtvc.npz"
                file_path = os.path.join(cell_dir, file_pattern)
                
                # 尝试加载染色体矩阵
                matrix = self._load_single_chromosome(file_path, chrom)
                
                if matrix is not None:
                    available_chroms.append(chrom)
                    chrom_matrices.append(matrix)
                    logger.info(f"  ✅ {chrom}: {matrix.shape}, 非零: {matrix.nnz:,}")
                    
                    # 内存管理：及时释放大矩阵的引用
                    if not use_memory_efficient:
                        pass  # 保持所有矩阵在内存中
                    else:
                        # 在内存高效模式下，暂时保持矩阵，后续会处理
                        pass
                        
                else:
                    # 创建空矩阵作为占位符
                    expected_bins = self.chrom_bins[chrom]
                    empty_matrix = csr_matrix((expected_bins, expected_bins))
                    available_chroms.append(chrom)
                    chrom_matrices.append(empty_matrix)
                    logger.warning(f"  ⚠️  {chrom}: 使用空矩阵代替")
                
                # 强制垃圾回收
                if use_memory_efficient:
                    gc.collect()
            
            # 检查是否有足够的染色体数据
            if len(available_chroms) < 10:  # 至少需要一半的染色体
                logger.error(f"可用染色体数量不足: {len(available_chroms)}/20")
                return False
            
            logger.info(f"找到 {len(available_chroms)} 个染色体数据")
            
            # 使用block_diag合并矩阵
            logger.info("开始矩阵合并...")
            merge_start = time.time()
            
            try:
                merged_matrix = block_diag(chrom_matrices, format='csr')
                merge_time = time.time() - merge_start
                
                logger.info(f"矩阵合并完成: {merged_matrix.shape}, 非零: {merged_matrix.nnz:,}, 耗时: {merge_time:.2f}秒")
                
                # 保存合并后的矩阵
                save_npz(output_path, merged_matrix)
                logger.info(f"矩阵保存完成: {output_path}")
                
                # 清理内存
                del merged_matrix
                del chrom_matrices
                gc.collect()
                
                elapsed_time = time.time() - start_time
                logger.info(f"细胞 {cell_name} 合并成功，总耗时: {elapsed_time:.2f}秒")
                
                return True
                
            except MemoryError as e:
                logger.error(f"内存不足，无法合并矩阵: {e}")
                # 尝试清理内存后重试
                del chrom_matrices
                gc.collect()
                return False
                
        except Exception as e:
            logger.error(f"合并细胞 {cell_name} 时出现异常: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return False
            
    def merge_stage_cells_robust(self, stage_dir: str, stage_name: str, 
                                output_dir: str, max_cells: Optional[int] = None,
                                use_memory_efficient: bool = True) -> Dict[str, str]:
        """
        鲁棒的stage细胞合并
        
        Parameters:
        -----------
        stage_dir : str
            stage目录路径
        stage_name : str
            stage名称
        output_dir : str
            输出目录
        max_cells : int, optional
            最大处理细胞数
        use_memory_efficient : bool
            是否使用内存高效模式
            
        Returns:
        --------
        Dict[str, str] : 成功合并的细胞文件映射
        """
        logger.info(f"开始鲁棒合并Stage {stage_name} 的细胞矩阵")
        
        # 创建输出目录
        stage_output_dir = os.path.join(output_dir, stage_name)
        os.makedirs(stage_output_dir, exist_ok=True)
        
        # 获取所有细胞目录
        cell_dirs = [d for d in os.listdir(stage_dir) 
                    if os.path.isdir(os.path.join(stage_dir, d)) and 
                    (d.startswith("Gas") or d.startswith("Org"))]
        
        if max_cells:
            cell_dirs = cell_dirs[:max_cells]
            
        logger.info(f"Stage {stage_name} 找到 {len(cell_dirs)} 个细胞目录")
        
        # 合并每个细胞
        merged_files = {}
        failed_cells = []
        
        for i, cell_name in enumerate(cell_dirs, 1):
            logger.info(f"处理细胞 {i}/{len(cell_dirs)}: {cell_name}")
            
            cell_dir = os.path.join(stage_dir, cell_name)
            output_path = os.path.join(stage_output_dir, f"{cell_name}_merged.npz")
            
            # 使用鲁棒合并方法
            success = self.merge_cell_chromosomes_robust(
                cell_dir, cell_name, output_path, use_memory_efficient
            )
            
            if success:
                merged_files[cell_name] = output_path
                logger.info(f"✅ 细胞 {cell_name} 合并成功")
            else:
                failed_cells.append(cell_name)
                logger.error(f"❌ 细胞 {cell_name} 合并失败")
                
            # 每处理几个细胞就强制垃圾回收
            if i % 3 == 0:
                gc.collect()
                
        logger.info(f"Stage {stage_name} 合并完成: {len(merged_files)} 成功, {len(failed_cells)} 失败")
        
        if failed_cells:
            logger.warning(f"失败的细胞: {failed_cells[:5]}{'...' if len(failed_cells) > 5 else ''}")
            
        return merged_files


def main():
    """测试主函数"""
    # 配置路径
    chrom_sizes_file = "/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes.txt"
    input_base_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage"
    output_base_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/robust_merged_matrices"
    
    # 创建鲁棒合并器
    merger = RobustHiCMatrixMerger(chrom_sizes_file, resolution=40000)
    
    # 测试较小的stage
    test_stages = ['E85', 'E95']  # 选择细胞数较少的stage进行测试
    
    for stage in test_stages:
        stage_dir = os.path.join(input_base_dir, stage)
        if os.path.exists(stage_dir):
            logger.info(f"开始测试Stage: {stage}")
            merged_files = merger.merge_stage_cells_robust(
                stage_dir, stage, output_base_dir, 
                max_cells=2,  # 先测试2个细胞
                use_memory_efficient=True
            )
            logger.info(f"Stage {stage} 测试完成: {len(merged_files)} 个细胞成功")
        else:
            logger.warning(f"Stage目录不存在: {stage_dir}")


if __name__ == "__main__":
    main()