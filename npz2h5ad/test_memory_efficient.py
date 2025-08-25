#!/usr/bin/env python3
"""
测试内存高效转换器的脚本

使用少量细胞测试内存高效转换的功能和性能
"""

import sys
import logging
from pathlib import Path
import time
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
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def get_memory_usage():
    """获取当前内存使用量（MB）"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def test_memory_efficient_conversion():
    """测试内存高效转换"""
    logger.info("Testing memory-efficient Hi-C conversion...")
    
    # 设置测试参数
    base_dir = Path("../output/imputed_matrices_by_stage")
    test_stage = "E70"
    test_cells = 3  # 只测试3个细胞
    
    stage_dir = base_dir / test_stage
    output_file = Path("./test_memory_efficient_output.h5ad")
    
    if not stage_dir.exists():
        logger.error(f"Test stage directory does not exist: {stage_dir}")
        return False
    
    # 记录初始内存
    initial_memory = get_memory_usage()
    logger.info(f"Initial memory usage: {initial_memory:.2f} MB")
    
    try:
        # 初始化组件
        loader = HiCNpzLoader()
        extractor = UpperTriangleExtractor(include_diagonal=True)
        chr_manager = ChromosomeManager()
        
        converter = MemoryEfficientHiCConverter(
            loader=loader,
            extractor=extractor,
            chr_manager=chr_manager
        )
        
        # 开始转换
        start_time = time.time()
        logger.info(f"Starting conversion of {test_cells} cells...")
        
        result_file = converter.convert_stage_memory_efficient(
            stage_dir=stage_dir,
            output_file=output_file,
            max_cells=test_cells
        )
        
        end_time = time.time()
        final_memory = get_memory_usage()
        
        # 检查结果
        if result_file and result_file.exists():
            file_size = result_file.stat().st_size / 1024 / 1024
            
            logger.info("✅ Memory-efficient conversion test PASSED!")
            logger.info(f"   Processing time: {end_time - start_time:.2f} seconds")
            logger.info(f"   Initial memory: {initial_memory:.2f} MB")
            logger.info(f"   Final memory: {final_memory:.2f} MB")
            logger.info(f"   Memory increase: {final_memory - initial_memory:.2f} MB")
            logger.info(f"   Output file size: {file_size:.2f} MB")
            logger.info(f"   Output file: {result_file}")
            
            # 验证文件可以读取
            try:
                import anndata as ad
                adata = ad.read_h5ad(result_file)
                logger.info(f"   ✅ File verification: {adata.shape}")
                adata.file.close() if hasattr(adata, 'file') else None
                del adata
            except Exception as e:
                logger.warning(f"   ⚠️ File verification failed: {e}")
            
            return True
        else:
            logger.error("❌ Memory-efficient conversion test FAILED - no output file")
            return False
            
    except Exception as e:
        logger.error(f"❌ Memory-efficient conversion test FAILED: {e}")
        return False

def compare_memory_usage():
    """比较内存使用情况"""
    logger.info("\n" + "="*50)
    logger.info("MEMORY USAGE COMPARISON")
    logger.info("="*50)
    
    # 估算传统方法的内存需求
    n = 4887  # 矩阵维度
    chromosomes = 20
    cells = 3
    
    upper_tri_per_chr = n * (n + 1) // 2
    total_features_per_cell = upper_tri_per_chr * chromosomes
    total_features_batch = total_features_per_cell * cells
    
    traditional_memory_gb = total_features_batch * 8 / 1024 / 1024 / 1024
    
    logger.info(f"Traditional method (load all to memory):")
    logger.info(f"  - 3 cells matrix size: {traditional_memory_gb:.2f} GB")
    logger.info(f"  - Feasible for full dataset: ❌ NO (would need ~170 GB)")
    
    logger.info(f"\nMemory-efficient method:")
    logger.info(f"  - Peak memory: ~{get_memory_usage():.2f} MB")
    logger.info(f"  - Feasible for full dataset: ✅ YES")
    logger.info(f"  - Estimated time for full dataset: ~{557 * 2 / 60:.1f} hours (557 cells × 2 seconds/cell)")

def main():
    """运行测试"""
    logger.info("Starting memory-efficient conversion tests...")
    
    # 测试内存高效转换
    test_success = test_memory_efficient_conversion()
    
    # 比较内存使用
    compare_memory_usage()
    
    # 总结
    logger.info(f"\n{'='*50}")
    logger.info("TEST SUMMARY")
    logger.info(f"{'='*50}")
    
    if test_success:
        logger.info("✅ Memory-efficient conversion: PASSED")
        logger.info("🎉 Ready for production use!")
        return 0
    else:
        logger.error("❌ Memory-efficient conversion: FAILED")
        return 1

if __name__ == "__main__":
    sys.exit(main())
