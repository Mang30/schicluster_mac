#!/usr/bin/env python3
"""
Hi-C NPZ 数据批量转换为 AnnData 格式的主程序

此脚本提供了完整的流水线，将组织在分期目录中的单细胞 Hi-C NPZ 数据
转换为 AnnData 格式，便于后续的单细胞基因组学分析。

使用方法:
    python convert_npz_to_h5ad.py --input_dir ../output/imputed_matrices_by_stage \
                                  --output_dir ./output \
                                  --stages E70,E75,E80 \
                                  --include_diagonal
"""

import argparse
import sys
from pathlib import Path
import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Optional
import warnings

# 导入我们的工具库
from hic_converter import (
    HiCNpzLoader, 
    UpperTriangleExtractor, 
    ChromosomeManager,
    FeatureGenerator, 
    AnnDataBuilder
)

# 设置日志
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('npz_to_h5ad_conversion.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Convert Hi-C NPZ data to AnnData format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # 转换所有阶段的数据
    python convert_npz_to_h5ad.py --input_dir ../output/imputed_matrices_by_stage --output_dir ./output
    
    # 只转换特定阶段
    python convert_npz_to_h5ad.py --input_dir ../output/imputed_matrices_by_stage --output_dir ./output --stages E70,E75
    
    # 不包含对角线元素
    python convert_npz_to_h5ad.py --input_dir ../output/imputed_matrices_by_stage --output_dir ./output --no_diagonal
        """
    )
    
    parser.add_argument(
        '--input_dir', 
        type=str, 
        default='../output/imputed_matrices_by_stage',
        help='Input directory containing stage folders (default: ../output/imputed_matrices_by_stage)'
    )
    
    parser.add_argument(
        '--output_dir', 
        type=str, 
        default='./output',
        help='Output directory for H5AD files (default: ./output)'
    )
    
    parser.add_argument(
        '--stages', 
        type=str, 
        default=None,
        help='Comma-separated list of stages to process (e.g., E70,E75,E80). If not specified, all stages will be processed.'
    )
    
    parser.add_argument(
        '--include_diagonal', 
        action='store_true',
        default=True,
        help='Include diagonal elements in upper triangular extraction (default: True)'
    )
    
    parser.add_argument(
        '--no_diagonal', 
        action='store_true',
        help='Exclude diagonal elements from upper triangular extraction'
    )
    
    parser.add_argument(
        '--file_pattern', 
        type=str, 
        default='*_{}_pad1_std1.0_rp0.5_sqrtvc.npz',
        help='File pattern for chromosome files (default: *_{}_pad1_std1.0_rp0.5_sqrtvc.npz)'
    )
    
    parser.add_argument(
        '--max_cells_per_stage', 
        type=int, 
        default=None,
        help='Maximum number of cells to process per stage (for testing)'
    )
    
    parser.add_argument(
        '--chromosomes', 
        type=str, 
        default=None,
        help='Comma-separated list of chromosomes to include (default: chr1-chr19,chrX)'
    )
    
    parser.add_argument(
        '--verbose', 
        action='store_true',
        help='Enable verbose logging'
    )
    
    return parser.parse_args()

def discover_stages(input_dir: Path) -> List[str]:
    """
    发现输入目录中的所有阶段
    
    Parameters:
    -----------
    input_dir : Path
        输入目录
        
    Returns:
    --------
    stages : list of str
        发现的阶段列表
    """
    stages = []
    for item in input_dir.iterdir():
        if item.is_dir() and item.name.startswith(('E', 'EX')):
            stages.append(item.name)
    
    stages.sort()
    logger.info(f"Discovered stages: {stages}")
    return stages

def discover_cells_in_stage(stage_dir: Path, max_cells: Optional[int] = None) -> List[Path]:
    """
    发现阶段目录中的所有细胞目录
    
    Parameters:
    -----------
    stage_dir : Path
        阶段目录
    max_cells : int, optional
        最大细胞数量限制
        
    Returns:
    --------
    cell_dirs : list of Path
        细胞目录列表
    """
    cell_dirs = []
    for item in stage_dir.iterdir():
        if item.is_dir():
            cell_dirs.append(item)
    
    cell_dirs.sort(key=lambda x: x.name)
    
    if max_cells is not None and len(cell_dirs) > max_cells:
        logger.info(f"Limiting to {max_cells} cells out of {len(cell_dirs)} found")
        cell_dirs = cell_dirs[:max_cells]
    
    logger.info(f"Found {len(cell_dirs)} cells in stage {stage_dir.name}")
    return cell_dirs

def process_single_stage(stage_name: str, 
                        stage_dir: Path,
                        output_dir: Path,
                        feature_generator: FeatureGenerator,
                        anndata_builder: AnnDataBuilder,
                        file_pattern: str,
                        max_cells: Optional[int] = None) -> Path:
    """
    处理单个发育阶段的数据
    
    Parameters:
    -----------
    stage_name : str
        阶段名称
    stage_dir : Path
        阶段目录
    output_dir : Path
        输出目录
    feature_generator : FeatureGenerator
        特征生成器
    anndata_builder : AnnDataBuilder
        AnnData构建器
    file_pattern : str
        文件模式
    max_cells : int, optional
        最大细胞数限制
        
    Returns:
    --------
    output_file : Path
        输出的H5AD文件路径
    """
    logger.info(f"Processing stage: {stage_name}")
    
    # 发现细胞目录
    cell_dirs = discover_cells_in_stage(stage_dir, max_cells)
    
    if not cell_dirs:
        logger.warning(f"No cells found in stage {stage_name}")
        return None
    
    # 处理每个细胞
    feature_matrices = []
    cell_ids = []
    feature_info = None
    
    for i, cell_dir in enumerate(cell_dirs):
        try:
            logger.info(f"Processing cell {i+1}/{len(cell_dirs)}: {cell_dir.name}")
            
            # 生成细胞特征
            cell_features, cell_feature_info = feature_generator.generate_cell_features(
                cell_dir, file_pattern
            )
            
            # 保存特征信息（所有细胞应该有相同的特征结构）
            if feature_info is None:
                feature_info = cell_feature_info
            elif len(feature_info) != len(cell_feature_info):
                logger.warning(f"Feature dimension mismatch for cell {cell_dir.name}: "
                             f"expected {len(feature_info)}, got {len(cell_feature_info)}")
                continue
            
            feature_matrices.append(cell_features)
            cell_ids.append(cell_dir.name)
            
        except Exception as e:
            logger.error(f"Failed to process cell {cell_dir.name}: {e}")
            continue
    
    if not feature_matrices:
        logger.error(f"No cells successfully processed for stage {stage_name}")
        return None
    
    # 组合成矩阵
    logger.info(f"Combining {len(feature_matrices)} cells into feature matrix")
    combined_matrix = np.vstack(feature_matrices)
    
    # 创建阶段信息
    stage_info = {cell_id: stage_name for cell_id in cell_ids}
    
    # 构建 AnnData 对象
    logger.info("Building AnnData object")
    adata = anndata_builder.create_anndata(
        feature_matrix=combined_matrix,
        cell_ids=cell_ids,
        feature_info=feature_info,
        stage_info=stage_info
    )
    
    # 保存文件
    output_file = output_dir / f"{stage_name}_hic_data.h5ad"
    logger.info(f"Saving AnnData to {output_file}")
    adata.write(output_file)
    
    # 输出统计信息
    logger.info(f"Stage {stage_name} summary:")
    logger.info(f"  Cells: {adata.n_obs}")
    logger.info(f"  Features: {adata.n_vars}")
    logger.info(f"  Chromosomes: {len(set(adata.var['chromosome']))}")
    logger.info(f"  Output file: {output_file}")
    
    return output_file

def main():
    """主函数"""
    args = parse_arguments()
    
    # 设置日志级别
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # 解析对角线参数
    include_diagonal = args.include_diagonal and not args.no_diagonal
    
    # 解析染色体列表
    if args.chromosomes:
        chromosomes = [chr.strip() for chr in args.chromosomes.split(',')]
    else:
        chromosomes = None
    
    logger.info("Starting Hi-C NPZ to AnnData conversion")
    logger.info(f"Input directory: {args.input_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Include diagonal: {include_diagonal}")
    logger.info(f"File pattern: {args.file_pattern}")
    
    # 验证输入目录
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        logger.error(f"Input directory does not exist: {input_dir}")
        sys.exit(1)
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 发现或解析阶段
    if args.stages:
        stages = [stage.strip() for stage in args.stages.split(',')]
        logger.info(f"Processing specified stages: {stages}")
    else:
        stages = discover_stages(input_dir)
        logger.info(f"Processing all discovered stages: {stages}")
    
    # 初始化组件
    loader = HiCNpzLoader()
    extractor = UpperTriangleExtractor(include_diagonal=include_diagonal)
    chr_manager = ChromosomeManager(chromosomes=chromosomes)
    feature_generator = FeatureGenerator(loader, extractor, chr_manager)
    anndata_builder = AnnDataBuilder()
    
    # 处理每个阶段
    output_files = []
    
    for stage_name in stages:
        stage_dir = input_dir / stage_name
        
        if not stage_dir.exists():
            logger.warning(f"Stage directory does not exist: {stage_dir}")
            continue
        
        try:
            output_file = process_single_stage(
                stage_name=stage_name,
                stage_dir=stage_dir,
                output_dir=output_dir,
                feature_generator=feature_generator,
                anndata_builder=anndata_builder,
                file_pattern=args.file_pattern,
                max_cells=args.max_cells_per_stage
            )
            
            if output_file:
                output_files.append(output_file)
                
        except Exception as e:
            logger.error(f"Failed to process stage {stage_name}: {e}")
            continue
    
    # 总结
    logger.info("Conversion completed!")
    logger.info(f"Successfully processed {len(output_files)} stages:")
    for output_file in output_files:
        logger.info(f"  {output_file}")
    
    if len(output_files) == 0:
        logger.error("No stages were successfully processed")
        sys.exit(1)

if __name__ == "__main__":
    main()
