#!/usr/bin/env python3
"""
JSON to H5AD 核心转换脚本
将contact decay profile从JSON格式转换为H5AD格式
"""

import argparse
import sys
import warnings
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

# 科学计算库
import numpy as np
import pandas as pd
import scanpy as sc

# 本地模块
from metadata_processor import MetadataProcessor
from data_validator import DataValidator
from utils import (
    load_config, setup_logging, ProgressTracker, MemoryMonitor,
    create_directory, backup_file, get_file_list, save_json,
    format_file_size, check_disk_space, get_system_info
)

# 抑制警告
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)


class JSON2H5ADConverter:
    """JSON到H5AD转换器"""
    
    def __init__(self, config_file: str = "config.yaml"):
        """
        初始化转换器
        
        Args:
            config_file: 配置文件路径
        """
        self.config = load_config(config_file)
        self.logger = None  # 将在process_stage中初始化
        self.metadata_processor = None
        self.data_validator = DataValidator()
        self.memory_monitor = MemoryMonitor()
        
        # 验证配置
        self._validate_config()
    
    def _validate_config(self) -> None:
        """验证配置文件"""
        from utils import ConfigValidator
        
        errors = ConfigValidator.validate_config(self.config)
        if errors:
            print("配置文件验证失败:")
            for error in errors:
                print(f"  - {error}")
            sys.exit(1)
    
    def process_stage(self, stage: str, debug_mode: bool = False) -> bool:
        """
        处理单个stage
        
        Args:
            stage: 发育阶段名称 (如 'E75')
            debug_mode: 是否为调试模式
            
        Returns:
            是否成功处理
        """
        # 初始化日志
        self.logger = setup_logging(self.config, stage)
        self.memory_monitor.logger = self.logger
        
        self.logger.info(f"开始处理 Stage {stage}")
        self.logger.info(f"系统信息: {get_system_info()}")
        
        try:
            # 检查输入目录
            stage_input_dir = Path(self.config['paths']['input_root']) / f"stage_{stage}"
            if not stage_input_dir.exists():
                self.logger.error(f"输入目录不存在: {stage_input_dir}")
                return False
            
            # 检查输出目录和磁盘空间
            output_dir = Path(self.config['paths']['output_dir'])
            create_directory(output_dir)
            
            if not check_disk_space(output_dir, required_gb=1.0):
                self.logger.error("磁盘空间不足")
                return False
            
            # 初始化元数据处理器
            metadata_file = self.config['paths']['metadata_file']
            self.metadata_processor = MetadataProcessor(metadata_file)
            
            # 获取JSON文件列表
            json_files = get_file_list(stage_input_dir, "*_decay_profile.json")
            if not json_files:
                self.logger.warning(f"Stage {stage} 没有找到JSON文件")
                return True
            
            # 调试模式限制文件数量
            if debug_mode:
                max_files = self.config.get('debug', {}).get('max_files_debug', 5)
                json_files = json_files[:max_files]
                self.logger.info(f"调试模式: 仅处理前 {len(json_files)} 个文件")
            
            self.logger.info(f"找到 {len(json_files)} 个JSON文件")
            
            # 处理JSON文件
            with self.memory_monitor.monitor_block(f"处理 Stage {stage} 的 {len(json_files)} 个文件"):
                stage_data = self._process_json_files(json_files, stage)
            
            if not stage_data:
                self.logger.error(f"Stage {stage} 没有有效数据")
                return False
            
            # 验证和对齐数据
            with self.memory_monitor.monitor_block("验证和对齐数据"):
                success = self._validate_and_align_data(stage_data, stage)
                if not success:
                    return False
            
            # 转换为H5AD
            with self.memory_monitor.monitor_block("转换为H5AD格式"):
                adata = self._convert_to_h5ad(stage_data, stage)
                if adata is None:
                    return False
            
            # 保存H5AD文件
            with self.memory_monitor.monitor_block("保存H5AD文件"):
                success = self._save_h5ad(adata, stage)
            
            self.memory_monitor.log_memory_usage("处理完成")
            return success
            
        except Exception as e:
            self.logger.error(f"处理 Stage {stage} 时发生错误: {e}", exc_info=True)
            return False
    
    def _process_json_files(self, json_files: List[Path], stage: str) -> List[Dict]:
        """
        处理JSON文件列表
        
        Args:
            json_files: JSON文件路径列表
            stage: stage名称
            
        Returns:
            处理后的数据列表
        """
        stage_data = []
        progress = ProgressTracker(
            len(json_files), 
            f"处理 Stage {stage} JSON文件",
            self.config.get('performance', {}).get('progress_update_interval', 10)
        )
        
        batch_size = self.config['processing']['batch_size']
        
        for i in range(0, len(json_files), batch_size):
            batch_files = json_files[i:i + batch_size]
            batch_data = self._process_batch(batch_files)
            stage_data.extend(batch_data)
            
            progress.update(len(batch_files))
            
            # 定期报告内存使用情况
            if i % (batch_size * 5) == 0 and i > 0:
                self.memory_monitor.log_memory_usage(f"处理了 {i + len(batch_files)} 个文件")
        
        progress.finish()
        self.logger.info(f"成功处理 {len(stage_data)} 个文件")
        
        return stage_data
    
    def _process_batch(self, batch_files: List[Path]) -> List[Dict]:
        """
        处理一批JSON文件
        
        Args:
            batch_files: 文件路径列表
            
        Returns:
            处理后的数据列表
        """
        batch_data = []
        
        for json_file in batch_files:
            try:
                # 验证和读取JSON文件
                is_valid, data = self.data_validator.validate_json_file(str(json_file))
                if not is_valid:
                    continue
                
                # 提取细胞ID
                cell_id = self.metadata_processor.extract_cell_id_from_filename(json_file.name)
                data['cell_id'] = cell_id
                
                # 获取元数据
                metadata = self.metadata_processor.get_cell_metadata(cell_id)
                data.update(metadata)
                
                batch_data.append(data)
                
            except Exception as e:
                self.logger.error(f"处理文件 {json_file} 时出错: {e}")
                continue
        
        return batch_data
    
    def _validate_and_align_data(self, stage_data: List[Dict], stage: str) -> bool:
        """
        验证和对齐stage数据
        
        Args:
            stage_data: stage数据列表
            stage: stage名称
            
        Returns:
            是否成功
        """
        self.data_validator.reset_validation_status()
        
        # 验证stage一致性
        _, stats = self.data_validator.validate_stage_consistency(stage_data, stage)
        
        self.logger.info(f"Stage {stage} 统计信息:")
        self.logger.info(f"  细胞数量: {stats['cell_count']}")
        self.logger.info(f"  距离长度范围: {stats['distances_length_range']}")
        self.logger.info(f"  最常见长度: {stats['distances_length_mode']}")
        self.logger.info(f"  分辨率: {stats['resolution_unique']}")
        self.logger.info(f"  总contacts范围: {stats['total_contacts_range']}")
        self.logger.info(f"  平均contacts: {stats['total_contacts_mean']:.1f}")
        
        # 检查验证结果
        summary = self.data_validator.get_validation_summary()
        if summary['has_errors']:
            max_errors = self.config.get('validation', {}).get('max_errors', 0)
            if summary['error_count'] > max_errors:
                self.logger.error(f"验证错误数 ({summary['error_count']}) 超过限制 ({max_errors})")
                self.data_validator.log_validation_summary()
                return False
        
        self.data_validator.log_validation_summary()
        return True
    
    def _convert_to_h5ad(self, stage_data: List[Dict], stage: str) -> Optional[sc.AnnData]:
        """
        将数据转换为AnnData格式
        
        Args:
            stage_data: stage数据列表
            stage: stage名称
            
        Returns:
            AnnData对象或None
        """
        if not stage_data:
            self.logger.error("没有数据可转换")
            return None
        
        try:
            # 对齐distances
            unified_distances, unified_distance_kb = self.data_validator.align_distances(stage_data)
            if len(unified_distances) == 0:
                self.logger.error("无法对齐distances")
                return None
            
            # 创建decay矩阵
            decay_matrix, cell_ids = self.data_validator.create_decay_matrix(stage_data, unified_distances)
            
            self.logger.info(f"创建了 {decay_matrix.shape[0]} × {decay_matrix.shape[1]} 的decay矩阵")
            
            # 创建obs (细胞元数据)
            obs_data = []
            for data in stage_data:
                cell_metadata = {}
                obs_fields = self.config['h5ad_structure']['obs_fields']
                
                for field, dtype in obs_fields.items():
                    value = data.get(field)
                    if value is not None and not pd.isna(value):
                        cell_metadata[field] = value
                    else:
                        # 使用默认值
                        if 'int' in dtype:
                            cell_metadata[field] = 0
                        elif 'float' in dtype:
                            cell_metadata[field] = np.nan
                        else:
                            cell_metadata[field] = self.config['defaults'].get(f"unknown_{field}", "Unknown")
                
                obs_data.append(cell_metadata)
            
            obs_df = pd.DataFrame(obs_data, index=cell_ids)
            
            # 创建var (特征元数据)
            var_data = {
                'distance_bins': unified_distances,
                'distance_kb': unified_distance_kb
            }
            var_df = pd.DataFrame(var_data, index=[f"bin_{i}" for i in range(len(unified_distances))])
            
            # 创建AnnData对象
            adata = sc.AnnData(
                X=decay_matrix,
                obs=obs_df,
                var=var_df,
                dtype=np.float32
            )
            
            # 添加uns信息
            adata.uns['processing_date'] = datetime.now().isoformat()
            adata.uns['stage'] = stage
            adata.uns['stage_summary'] = self.metadata_processor.get_stage_summary().get(stage.upper(), {})
            adata.uns['celltype_counts'] = self.metadata_processor.get_celltype_counts(stage)
            adata.uns['processing_parameters'] = {
                'alignment_strategy': self.config['processing']['alignment_strategy'],
                'data_dtype': str(decay_matrix.dtype),
                'unified_distance_length': len(unified_distances)
            }
            adata.uns['validation_summary'] = self.data_validator.get_validation_summary()
            
            self.logger.info(f"创建AnnData对象: {adata.n_obs} 细胞 × {adata.n_vars} 特征")
            
            return adata
            
        except Exception as e:
            self.logger.error(f"转换为AnnData时出错: {e}", exc_info=True)
            return None
    
    def _save_h5ad(self, adata: sc.AnnData, stage: str) -> bool:
        """
        保存H5AD文件
        
        Args:
            adata: AnnData对象
            stage: stage名称
            
        Returns:
            是否成功保存
        """
        try:
            output_dir = Path(self.config['paths']['output_dir'])
            filename_format = self.config['output']['h5ad_filename_format']
            filename = filename_format.format(stage=stage)
            output_file = output_dir / filename
            
            # 检查是否需要备份
            if output_file.exists():
                if not self.config['output']['overwrite_existing']:
                    self.logger.error(f"输出文件已存在且不允许覆盖: {output_file}")
                    return False
                
                if self.config['output']['create_backup']:
                    backup_path = backup_file(str(output_file))
                    if backup_path:
                        self.logger.info(f"已创建备份文件: {backup_path}")
            
            # 保存H5AD文件
            self.logger.info(f"保存H5AD文件: {output_file}")
            
            if self.config['output']['compress_output']:
                adata.write_h5ad(
                    output_file,
                    compression='gzip',
                    compression_opts=self.config['output']['compression_level']
                )
            else:
                adata.write_h5ad(output_file)
            
            # 检查文件大小
            file_size = output_file.stat().st_size
            self.logger.info(f"H5AD文件已保存: {output_file} ({format_file_size(file_size)})")
            
            # 保存处理摘要
            self._save_processing_summary(adata, stage)
            
            return True
            
        except Exception as e:
            self.logger.error(f"保存H5AD文件时出错: {e}", exc_info=True)
            return False
    
    def _save_processing_summary(self, adata: sc.AnnData, stage: str) -> None:
        """
        保存处理摘要
        
        Args:
            adata: AnnData对象
            stage: stage名称
        """
        try:
            output_dir = Path(self.config['paths']['output_dir'])
            summary_file = output_dir / f"stage_{stage}_summary.json"
            
            summary = {
                'stage': stage,
                'processing_date': datetime.now().isoformat(),
                'cell_count': adata.n_obs,
                'feature_count': adata.n_vars,
                'celltype_distribution': adata.obs['celltype'].value_counts().to_dict(),
                'data_shape': adata.shape,
                'memory_usage_mb': adata.X.nbytes / (1024 * 1024),
                'validation_summary': adata.uns.get('validation_summary', {}),
                'processing_parameters': adata.uns.get('processing_parameters', {}),
                'system_info': get_system_info()
            }
            
            save_json(summary, str(summary_file))
            self.logger.info(f"处理摘要已保存: {summary_file}")
            
        except Exception as e:
            self.logger.warning(f"保存处理摘要时出错: {e}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="JSON to H5AD 转换器")
    parser.add_argument('--stage', type=str, required=True,
                       help='发育阶段 (E70, E75, E80, E85, E95, EX05, EX15)')
    parser.add_argument('--config', type=str, default='config.yaml',
                       help='配置文件路径')
    parser.add_argument('--debug', action='store_true',
                       help='调试模式 (处理少量文件)')
    parser.add_argument('--metadata', type=str,
                       help='元数据文件路径 (覆盖配置文件中的设置)')
    
    args = parser.parse_args()
    
    # 验证stage参数
    valid_stages = ['E70', 'E75', 'E80', 'E85', 'E95', 'EX05', 'EX15']
    if args.stage not in valid_stages:
        print(f"错误: 无效的stage '{args.stage}'. 有效值: {valid_stages}")
        sys.exit(1)
    
    try:
        # 创建转换器
        converter = JSON2H5ADConverter(args.config)
        
        # 如果指定了元数据文件，更新配置
        if args.metadata:
            converter.config['paths']['metadata_file'] = args.metadata
        
        # 处理stage
        success = converter.process_stage(args.stage, args.debug)
        
        if success:
            print(f"Stage {args.stage} 处理成功!")
            sys.exit(0)
        else:
            print(f"Stage {args.stage} 处理失败!")
            sys.exit(1)
            
    except Exception as e:
        print(f"发生错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()