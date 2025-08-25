#!/usr/bin/env python3
"""
内存高效转换的主程序

使用内存高效的方法转换 Hi-C NPZ 数据到 AnnData 格式
"""

import argparse
import sys
from pathlib import Path
import logging
import psutil
import os

# 导入工具库
from hic_converter import (
    HiCNpzLoader, 
    UpperTriangleExtractor, 
    ChromosomeManager
)
from memory_efficient_converter import MemoryEfficientHiCConverter

# 设置日志
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('memory_efficient_conversion.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def get_memory_usage():
    """获取当前内存使用量"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Memory-efficient Hi-C NPZ to AnnData conversion',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # 转换单个阶段（测试）
    python convert_memory_efficient.py --stage E70 --max_cells 10
    
    # 转换单个阶段（完整）
    python convert_memory_efficient.py --stage E70
    
    # 转换所有阶段
    python convert_memory_efficient.py --all_stages
        """
    )
    
    parser.add_argument(
        '--input_dir', 
        type=str, 
        default='../output/imputed_matrices_by_stage',
        help='Input directory containing stage folders'
    )
    
    parser.add_argument(
        '--output_dir', 
        type=str, 
        default='./output',
        help='Output directory for H5AD files'
    )
    
    parser.add_argument(
        '--stage', 
        type=str, 
        help='Single stage to process (e.g., E70)'
    )
    
    parser.add_argument(
        '--all_stages', 
        action='store_true',
        help='Process all discovered stages'
    )
    
    parser.add_argument(
        '--max_cells', 
        type=int, 
        default=None,
        help='Maximum number of cells to process (for testing)'
    )
    
    parser.add_argument(
        '--include_diagonal', 
        action='store_true',
        default=True,
        help='Include diagonal elements (default: True)'
    )
    
    parser.add_argument(
        '--no_diagonal', 
        action='store_true',
        help='Exclude diagonal elements'
    )
    
    parser.add_argument(
        '--chromosomes', 
        type=str, 
        default=None,
        help='Comma-separated list of chromosomes (e.g., 1,2,3,X)'
    )
    
    parser.add_argument(
        '--verbose', 
        action='store_true',
        help='Enable verbose logging'
    )
    
    return parser.parse_args()

def discover_stages(input_dir: Path) -> list:
    """发现所有阶段"""
    stages = []
    for item in input_dir.iterdir():
        if item.is_dir() and item.name.startswith(('E', 'EX')):
            stages.append(item.name)
    stages.sort()
    return stages

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
    
    logger.info("Starting memory-efficient Hi-C NPZ to AnnData conversion")
    logger.info(f"Input directory: {args.input_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Include diagonal: {include_diagonal}")
    logger.info(f"Initial memory usage: {get_memory_usage():.2f} MB")
    
    # 验证输入目录
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        logger.error(f"Input directory does not exist: {input_dir}")
        sys.exit(1)
    
    # 创建输出目录
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 确定要处理的阶段
    if args.stage:
        stages = [args.stage]
        logger.info(f"Processing single stage: {args.stage}")
    elif args.all_stages:
        stages = discover_stages(input_dir)
        logger.info(f"Processing all stages: {stages}")
    else:
        logger.error("Must specify either --stage or --all_stages")
        sys.exit(1)
    
    # 初始化组件
    loader = HiCNpzLoader()
    extractor = UpperTriangleExtractor(include_diagonal=include_diagonal)
    chr_manager = ChromosomeManager(chromosomes=chromosomes)
    
    converter = MemoryEfficientHiCConverter(
        loader=loader,
        extractor=extractor,
        chr_manager=chr_manager
    )
    
    # 处理每个阶段
    success_count = 0
    
    for stage_name in stages:
        stage_dir = input_dir / stage_name
        
        if not stage_dir.exists():
            logger.warning(f"Stage directory does not exist: {stage_dir}")
            continue
        
        output_file = output_dir / f"{stage_name}_hic_memory_efficient.h5ad"
        
        try:
            logger.info(f"\n{'='*60}")
            logger.info(f"Processing stage: {stage_name}")
            logger.info(f"Memory before stage: {get_memory_usage():.2f} MB")
            logger.info(f"{'='*60}")
            
            # 内存高效转换
            result_file = converter.convert_stage_memory_efficient(
                stage_dir=stage_dir,
                output_file=output_file,
                max_cells=args.max_cells
            )
            
            if result_file:
                success_count += 1
                file_size = result_file.stat().st_size / 1024 / 1024 / 1024
                logger.info(f"✅ Stage {stage_name} completed successfully")
                logger.info(f"   Output file: {result_file}")
                logger.info(f"   File size: {file_size:.2f} GB")
                logger.info(f"   Memory after stage: {get_memory_usage():.2f} MB")
            
        except Exception as e:
            logger.error(f"❌ Failed to process stage {stage_name}: {e}")
            continue
    
    # 总结
    logger.info(f"\n{'='*60}")
    logger.info("CONVERSION COMPLETED")
    logger.info(f"{'='*60}")
    logger.info(f"Successfully processed: {success_count}/{len(stages)} stages")
    logger.info(f"Final memory usage: {get_memory_usage():.2f} MB")
    
    if success_count == 0:
        logger.error("No stages were successfully processed")
        sys.exit(1)
    else:
        logger.info("🎉 Memory-efficient conversion completed successfully!")

if __name__ == "__main__":
    main()
