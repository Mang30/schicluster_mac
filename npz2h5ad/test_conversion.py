#!/usr/bin/env python3
"""
Hi-C NPZ 到 AnnData 转换工具的测试脚本

此脚本用于测试转换流程，使用小规模数据验证功能是否正常工作。
"""

import sys
import logging
from pathlib import Path
import numpy as np

# 导入我们的工具库
from hic_converter import (
    HiCNpzLoader, 
    UpperTriangleExtractor, 
    ChromosomeManager,
    FeatureGenerator, 
    AnnDataBuilder
)

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_single_cell():
    """测试单个细胞的处理流程"""
    logger.info("Testing single cell processing...")
    
    # 设置路径
    base_dir = Path("../output/imputed_matrices_by_stage")
    test_stage = "E70"
    test_cell = "GasdE701001"
    
    cell_dir = base_dir / test_stage / test_cell
    
    if not cell_dir.exists():
        logger.error(f"Test cell directory does not exist: {cell_dir}")
        return False
    
    # 初始化组件
    loader = HiCNpzLoader()
    extractor = UpperTriangleExtractor(include_diagonal=True)
    chr_manager = ChromosomeManager()
    feature_generator = FeatureGenerator(loader, extractor, chr_manager)
    
    try:
        # 测试特征生成，使用正确的文件模式
        file_pattern = "*_{}_pad1_std1.0_rp0.5_sqrtvc.npz"
        logger.info(f"Processing test cell: {test_cell}")
        features, feature_info = feature_generator.generate_cell_features(cell_dir, file_pattern)
        
        logger.info(f"Generated features: shape {features.shape}")
        logger.info(f"Feature count: {len(feature_info)}")
        logger.info(f"First 5 features: {feature_info[:5]}")
        logger.info(f"Feature value range: [{np.min(features):.6f}, {np.max(features):.6f}]")
        
        return True
        
    except Exception as e:
        logger.error(f"Failed to process test cell: {e}")
        return False

def test_small_batch():
    """测试小批量细胞处理"""
    logger.info("Testing small batch processing...")
    
    # 设置路径
    base_dir = Path("../output/imputed_matrices_by_stage")
    test_stage = "E70"
    stage_dir = base_dir / test_stage
    
    if not stage_dir.exists():
        logger.error(f"Test stage directory does not exist: {stage_dir}")
        return False
    
    # 获取前3个细胞进行测试
    cell_dirs = [d for d in stage_dir.iterdir() if d.is_dir()][:3]
    
    if len(cell_dirs) == 0:
        logger.error("No cell directories found for testing")
        return False
    
    logger.info(f"Testing with {len(cell_dirs)} cells: {[d.name for d in cell_dirs]}")
    
    # 初始化组件
    loader = HiCNpzLoader()
    extractor = UpperTriangleExtractor(include_diagonal=True)
    chr_manager = ChromosomeManager()
    feature_generator = FeatureGenerator(loader, extractor, chr_manager)
    anndata_builder = AnnDataBuilder()
    
    try:
        # 处理每个细胞，使用正确的文件模式
        feature_matrices = []
        cell_ids = []
        feature_info = None
        file_pattern = "*_{}_pad1_std1.0_rp0.5_sqrtvc.npz"
        
        for cell_dir in cell_dirs:
            logger.info(f"Processing cell: {cell_dir.name}")
            
            cell_features, cell_feature_info = feature_generator.generate_cell_features(cell_dir, file_pattern)
            
            if feature_info is None:
                feature_info = cell_feature_info
            
            feature_matrices.append(cell_features)
            cell_ids.append(cell_dir.name)
        
        # 组合矩阵
        combined_matrix = np.vstack(feature_matrices)
        logger.info(f"Combined matrix shape: {combined_matrix.shape}")
        
        # 创建 AnnData 对象
        stage_info = {cell_id: test_stage for cell_id in cell_ids}
        adata = anndata_builder.create_anndata(
            feature_matrix=combined_matrix,
            cell_ids=cell_ids,
            feature_info=feature_info,
            stage_info=stage_info
        )
        
        logger.info(f"Created AnnData object: {adata}")
        logger.info(f"Shape: {adata.shape}")
        logger.info(f"Observations: {adata.obs.head()}")
        logger.info(f"Variables info: {adata.var.head()}")
        
        # 保存测试文件
        test_output = Path("./test_output.h5ad")
        adata.write(test_output)
        logger.info(f"Test output saved to: {test_output}")
        
        return True
        
    except Exception as e:
        logger.error(f"Failed batch processing test: {e}")
        return False

def test_matrix_loading():
    """测试矩阵加载功能"""
    logger.info("Testing matrix loading...")
    
    # 测试加载样例数据
    example_files = [
        "../data/example_data/GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz",
        "../data/example_data/Astro_Endo_ODC_OPC_chr1_upper_tri.npz"
    ]
    
    loader = HiCNpzLoader()
    extractor = UpperTriangleExtractor(include_diagonal=True)
    
    for file_path in example_files:
        file_path = Path(file_path)
        if not file_path.exists():
            logger.warning(f"Example file does not exist: {file_path}")
            continue
        
        try:
            logger.info(f"Loading: {file_path.name}")
            matrix = loader.load_sparse_matrix(file_path)
            logger.info(f"Matrix shape: {matrix.shape}")
            logger.info(f"Matrix format: {type(matrix)}")
            logger.info(f"Non-zero elements: {matrix.nnz}")
            
            if matrix.shape[0] == matrix.shape[1]:  # 方形矩阵
                vector = extractor.extract_upper_triangle(matrix)
                logger.info(f"Upper triangle vector length: {len(vector)}")
                logger.info(f"Value range: [{np.min(vector):.6f}, {np.max(vector):.6f}]")
            else:
                logger.info("Non-square matrix, skipping upper triangle extraction")
            
        except Exception as e:
            logger.error(f"Failed to load {file_path}: {e}")
    
    return True

def main():
    """运行所有测试"""
    logger.info("Starting Hi-C conversion tool tests...")
    
    tests = [
        ("Matrix Loading", test_matrix_loading),
        ("Single Cell Processing", test_single_cell),
        ("Small Batch Processing", test_small_batch)
    ]
    
    results = {}
    
    for test_name, test_func in tests:
        logger.info(f"\n{'='*50}")
        logger.info(f"Running test: {test_name}")
        logger.info(f"{'='*50}")
        
        try:
            success = test_func()
            results[test_name] = success
            status = "PASSED" if success else "FAILED"
            logger.info(f"Test {test_name}: {status}")
        except Exception as e:
            logger.error(f"Test {test_name} crashed: {e}")
            results[test_name] = False
    
    # 总结测试结果
    logger.info(f"\n{'='*50}")
    logger.info("TEST SUMMARY")
    logger.info(f"{'='*50}")
    
    for test_name, success in results.items():
        status = "PASSED" if success else "FAILED"
        logger.info(f"{test_name:30s}: {status}")
    
    total_tests = len(results)
    passed_tests = sum(results.values())
    
    logger.info(f"\nTotal tests: {total_tests}")
    logger.info(f"Passed: {passed_tests}")
    logger.info(f"Failed: {total_tests - passed_tests}")
    
    if passed_tests == total_tests:
        logger.info("All tests PASSED! 🎉")
        return 0
    else:
        logger.error("Some tests FAILED! 😞")
        return 1

if __name__ == "__main__":
    sys.exit(main())
