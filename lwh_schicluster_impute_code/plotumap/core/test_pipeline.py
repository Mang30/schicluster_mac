#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hi-C UMAP可视化测试管道
基于E70, E80, EX05三个完整stage进行功能验证
"""

import os
import sys
import logging
import time
from pathlib import Path

# 添加当前目录到路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from matrix_merger_robust import RobustHiCMatrixMerger
from hic_anndata_builder import HiCAnnDataBuilder
from hic_umap_visualizer import HiCUMAPVisualizer

# 设置日志
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('hic_umap_pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class HiCUMAPPipeline:
    """Hi-C UMAP可视化完整管道"""
    
    def __init__(self, config: dict):
        """
        初始化管道
        
        Parameters:
        -----------
        config : dict
            配置参数字典
        """
        self.config = config
        self.validate_config()
        
        # 初始化各个模块
        self.merger = RobustHiCMatrixMerger(
            config['chrom_sizes_file'], 
            config['resolution']
        )
        
        self.builder = HiCAnnDataBuilder(
            config['chrom_sizes_file'], 
            config['resolution']
        )
        
        self.visualizer = HiCUMAPVisualizer(
            config['color_mapping_file']
        )
        
        logger.info("Hi-C UMAP Pipeline 初始化完成")
        
    def validate_config(self):
        """验证配置参数"""
        required_keys = [
            'chrom_sizes_file', 'input_base_dir', 'output_base_dir',
            'color_mapping_file', 'stage_files_mapping', 'resolution'
        ]
        
        for key in required_keys:
            if key not in self.config:
                raise ValueError(f"配置缺少必需参数: {key}")
                
        # 检查必要文件是否存在
        essential_files = [
            'chrom_sizes_file', 'color_mapping_file', 'stage_files_mapping'
        ]
        
        for file_key in essential_files:
            file_path = self.config[file_key]
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"必要文件不存在: {file_path}")
                
        logger.info("配置验证通过")
        
    def create_output_directories(self):
        """创建所有需要的输出目录"""
        output_dirs = [
            self.config['merged_matrices_dir'],
            self.config['h5ad_files_dir'],
            self.config['umap_plots_dir']
        ]
        
        for dir_path in output_dirs:
            os.makedirs(dir_path, exist_ok=True)
            logger.info(f"创建输出目录: {dir_path}")
            
    def step1_merge_matrices(self, stages: list, max_cells_per_stage: int = None) -> dict:
        """
        步骤1: 合并每个细胞的多染色体矩阵
        
        Parameters:
        -----------
        stages : list
            要处理的stage列表
        max_cells_per_stage : int, optional
            每个stage最大处理细胞数（用于测试）
            
        Returns:
        --------
        dict : 各stage的合并结果
        """
        logger.info("=== 步骤1: 开始矩阵合并 ===")
        merge_results = {}
        
        for stage in stages:
            stage_dir = os.path.join(self.config['input_base_dir'], stage)
            
            if not os.path.exists(stage_dir):
                logger.warning(f"Stage目录不存在: {stage_dir}")
                merge_results[stage] = {}
                continue
                
            logger.info(f"处理Stage: {stage}")
            start_time = time.time()
            
            merged_files = self.merger.merge_stage_cells_robust(
                stage_dir, 
                stage, 
                self.config['merged_matrices_dir'],
                max_cells_per_stage,
                use_memory_efficient=True
            )
            
            elapsed_time = time.time() - start_time
            merge_results[stage] = merged_files
            
            logger.info(f"Stage {stage} 合并完成: {len(merged_files)} 个细胞, "
                       f"耗时: {elapsed_time:.2f}秒")
                       
        total_merged = sum(len(files) for files in merge_results.values())
        logger.info(f"=== 步骤1完成: 总共合并 {total_merged} 个细胞矩阵 ===")
        
        return merge_results
        
    def step2_build_anndata(self, merge_results: dict) -> dict:
        """
        步骤2: 构建AnnData对象
        
        Parameters:
        -----------
        merge_results : dict
            矩阵合并结果
            
        Returns:
        --------
        dict : 各stage的AnnData构建结果
        """
        logger.info("=== 步骤2: 开始构建AnnData ===")
        
        # 加载metadata
        metadata_df = self.builder.load_cell_metadata(
            self.config['stage_files_mapping']
        )
        
        build_results = {}
        
        for stage, merged_files in merge_results.items():
            if not merged_files:
                logger.warning(f"Stage {stage} 没有合并后的矩阵文件")
                build_results[stage] = False
                continue
                
            logger.info(f"构建Stage {stage} 的AnnData")
            start_time = time.time()
            
            output_path = os.path.join(
                self.config['h5ad_files_dir'], 
                f"{stage}_hic.h5ad"
            )
            
            success = self.builder.build_stage_anndata(
                merged_files, 
                stage, 
                metadata_df, 
                output_path,
                normalization_method='log_cpm'
            )
            
            elapsed_time = time.time() - start_time
            build_results[stage] = success
            
            logger.info(f"Stage {stage} AnnData构建{'成功' if success else '失败'}, "
                       f"耗时: {elapsed_time:.2f}秒")
                       
        successful_builds = sum(build_results.values())
        logger.info(f"=== 步骤2完成: {successful_builds}/{len(build_results)} 个stage成功 ===")
        
        return build_results
        
    def step3_generate_umaps(self, build_results: dict) -> dict:
        """
        步骤3: 生成UMAP可视化
        
        Parameters:
        -----------
        build_results : dict
            AnnData构建结果
            
        Returns:
        --------
        dict : UMAP生成结果
        """
        logger.info("=== 步骤3: 开始UMAP可视化 ===")
        
        # 获取成功构建的stage列表
        successful_stages = [stage for stage, success in build_results.items() if success]
        
        if not successful_stages:
            logger.error("没有成功构建的AnnData，无法进行UMAP可视化")
            return {}
            
        logger.info(f"将为 {len(successful_stages)} 个stage生成UMAP: {successful_stages}")
        
        start_time = time.time()
        umap_results = self.visualizer.process_all_stages(
            self.config['h5ad_files_dir'],
            self.config['umap_plots_dir'],
            successful_stages
        )
        
        elapsed_time = time.time() - start_time
        successful_umaps = sum(umap_results.values())
        
        logger.info(f"=== 步骤3完成: {successful_umaps}/{len(umap_results)} 个stage成功, "
                   f"耗时: {elapsed_time:.2f}秒 ===")
                   
        return umap_results
        
    def run_full_pipeline(self, stages: list, max_cells_per_stage: int = None) -> dict:
        """
        运行完整的处理管道
        
        Parameters:
        -----------
        stages : list
            要处理的stage列表
        max_cells_per_stage : int, optional
            每个stage最大处理细胞数（用于测试）
            
        Returns:
        --------
        dict : 完整的处理结果
        """
        logger.info("🚀 开始运行Hi-C UMAP可视化完整管道")
        pipeline_start_time = time.time()
        
        # 创建输出目录
        self.create_output_directories()
        
        try:
            # 步骤1: 矩阵合并
            merge_results = self.step1_merge_matrices(stages, max_cells_per_stage)
            
            # 步骤2: AnnData构建
            build_results = self.step2_build_anndata(merge_results)
            
            # 步骤3: UMAP可视化
            umap_results = self.step3_generate_umaps(build_results)
            
            # 整理结果
            final_results = {
                'merge_results': merge_results,
                'build_results': build_results, 
                'umap_results': umap_results,
                'summary': {
                    'total_stages': len(stages),
                    'merged_cells': sum(len(files) for files in merge_results.values()),
                    'successful_anndata': sum(build_results.values()),
                    'successful_umaps': sum(umap_results.values())
                }
            }
            
            pipeline_elapsed_time = time.time() - pipeline_start_time
            
            logger.info("🎉 完整管道执行完成!")
            logger.info(f"📊 执行总结:")
            logger.info(f"   - 处理stage数: {final_results['summary']['total_stages']}")
            logger.info(f"   - 合并细胞数: {final_results['summary']['merged_cells']}")
            logger.info(f"   - 成功AnnData: {final_results['summary']['successful_anndata']}")
            logger.info(f"   - 成功UMAP: {final_results['summary']['successful_umaps']}")
            logger.info(f"   - 总耗时: {pipeline_elapsed_time:.2f}秒")
            
            return final_results
            
        except Exception as e:
            logger.error(f"管道执行失败: {e}")
            raise


def create_test_config() -> dict:
    """创建测试配置"""
    try:
        # 尝试使用配置文件
        from .config import CONFIG
        config = {
            # 输入配置
            'chrom_sizes_file': CONFIG.CHROM_SIZES,
            'input_base_dir': CONFIG.IMPUTED_MATRICES,
            'stage_files_mapping': CONFIG.STAGE_FILES_MAPPING,
            'color_mapping_file': CONFIG.COLOR_MAPPING,
            
            # 输出配置
            'output_base_dir': CONFIG.OUTPUT_DIR,
            'merged_matrices_dir': CONFIG.MERGED_MATRICES,
            'h5ad_files_dir': CONFIG.H5AD_OUTPUT,
            'umap_plots_dir': f"{CONFIG.OUTPUT_DIR}/umap_plots",
            
            # 处理参数
            'resolution': 40000  # 40kb (已通过实际数据检测确认)
        }
    except ImportError:
        # 回退到硬编码路径
        base_path = "/Volumes/SumSung500/CSU/0_HiRES"
        config = {
            # 输入配置
            'chrom_sizes_file': f"{base_path}/mm10_chrom_sizes.txt",
            'input_base_dir': f"{base_path}/output/imputed_matrices_by_stage",
            'stage_files_mapping': f"{base_path}/lwh_schicluster_impute_code/plotumap/config/stage_files_mapping.csv",
            'color_mapping_file': f"{base_path}/lwh_schicluster_impute_code/plotumap/config/color_mapping.json",
            
            # 输出配置
            'output_base_dir': f"{base_path}/output",
            'merged_matrices_dir': f"{base_path}/output/merged_matrices",
            'h5ad_files_dir': f"{base_path}/output/hic_h5ad_files",
            'umap_plots_dir': f"{base_path}/output/umap_plots",
            
            # 处理参数
            'resolution': 40000  # 40kb (已通过实际数据检测确认)
        }
    
    return config


def main():
    """主函数 - 运行测试管道"""
    logger.info("🧪 开始Hi-C UMAP可视化测试")
    
    # 创建配置
    config = create_test_config()
    
    # 测试阶段（只使用完整的stage）
    test_stages = ['E70', 'E80', 'EX05']
    
    # 测试参数（限制细胞数以便快速测试）
    max_cells_per_stage = 5  # 每个stage最多处理5个细胞
    
    try:
        # 创建管道
        pipeline = HiCUMAPPipeline(config)
        
        # 运行完整管道
        results = pipeline.run_full_pipeline(test_stages, max_cells_per_stage)
        
        logger.info("✅ 测试管道执行成功!")
        
        # 输出结果文件位置
        logger.info("📁 输出文件位置:")
        logger.info(f"   - 合并矩阵: {config['merged_matrices_dir']}")
        logger.info(f"   - AnnData文件: {config['h5ad_files_dir']}")
        logger.info(f"   - UMAP图表: {config['umap_plots_dir']}")
        
        return results
        
    except Exception as e:
        logger.error(f"❌ 测试管道执行失败: {e}")
        raise


if __name__ == "__main__":
    results = main()