#!/usr/bin/env python3
"""
数据验证和质量控制模块
用于验证JSON文件格式、数据一致性和质量控制
"""

import json
import logging
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import pandas as pd
from collections import Counter


class DataValidator:
    """数据验证器"""
    
    def __init__(self):
        """初始化数据验证器"""
        self.logger = logging.getLogger(__name__)
        self.validation_errors = []
        self.validation_warnings = []
    
    def reset_validation_status(self) -> None:
        """重置验证状态"""
        self.validation_errors.clear()
        self.validation_warnings.clear()
    
    def validate_json_file(self, json_file: str) -> Tuple[bool, Dict]:
        """
        验证单个JSON文件
        
        Args:
            json_file: JSON文件路径
            
        Returns:
            (是否有效, 数据字典或None)
        """
        json_path = Path(json_file)
        
        # 检查文件是否存在
        if not json_path.exists():
            error_msg = f"文件不存在: {json_file}"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            return False, {}
        
        # 检查文件名格式
        if not self._validate_filename(json_path.name):
            error_msg = f"文件名格式不正确: {json_path.name}"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            return False, {}
        
        try:
            # 读取JSON文件
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # 验证数据结构
            is_valid = self._validate_json_structure(data, json_file)
            
            if is_valid:
                # 验证数据内容
                is_valid = self._validate_json_content(data, json_file)
            
            return is_valid, data if is_valid else {}
            
        except json.JSONDecodeError as e:
            error_msg = f"JSON格式错误 {json_file}: {e}"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            return False, {}
        
        except Exception as e:
            error_msg = f"读取文件错误 {json_file}: {e}"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            return False, {}
    
    def _validate_filename(self, filename: str) -> bool:
        """
        验证文件名格式
        
        Args:
            filename: 文件名
            
        Returns:
            是否有效
        """
        # 跳过分析摘要文件
        if filename == 'analysis_summary.json':
            return False
        
        # 检查是否为decay profile文件
        if not filename.endswith('_decay_profile.json'):
            return False
        
        # 检查细胞ID格式 (基本格式检查)
        cell_id = filename.replace('_decay_profile.json', '')
        if len(cell_id) < 5:  # 最小长度检查
            return False
        
        return True
    
    def _validate_json_structure(self, data: Dict, filename: str) -> bool:
        """
        验证JSON数据结构
        
        Args:
            data: JSON数据
            filename: 文件名
            
        Returns:
            是否有效
        """
        required_fields = ['distances', 'decay_values', 'distance_kb', 
                          'total_contacts', 'resolution', 'power_law_slope']
        
        missing_fields = []
        for field in required_fields:
            if field not in data:
                missing_fields.append(field)
        
        if missing_fields:
            error_msg = f"文件 {filename} 缺少必需字段: {missing_fields}"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            return False
        
        return True
    
    def _validate_json_content(self, data: Dict, filename: str) -> bool:
        """
        验证JSON数据内容
        
        Args:
            data: JSON数据
            filename: 文件名
            
        Returns:
            是否有效
        """
        is_valid = True
        
        # 获取数据数组
        distances = data.get('distances', [])
        decay_values = data.get('decay_values', [])
        distance_kb = data.get('distance_kb', [])
        
        # 验证数组长度一致性
        lengths = [len(distances), len(decay_values), len(distance_kb)]
        if len(set(lengths)) > 1:
            error_msg = f"文件 {filename} 数据长度不一致: distances({lengths[0]}), decay_values({lengths[1]}), distance_kb({lengths[2]})"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            is_valid = False
        
        # 验证数据不为空
        if lengths[0] == 0:
            error_msg = f"文件 {filename} 数据为空"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            is_valid = False
            return is_valid
        
        # 验证distances是递增序列或多段递增序列
        if not self._is_increasing(distances):
            error_msg = f"文件 {filename} distances包含无效的非递增序列"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            is_valid = False
        
        # 验证decay_values为非负数
        if not self._all_non_negative(decay_values):
            error_msg = f"文件 {filename} decay_values包含负值"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            is_valid = False
        
        # 验证distance_kb为正数
        if not self._all_positive(distance_kb):
            error_msg = f"文件 {filename} distance_kb包含非正值"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            is_valid = False
        
        # 验证数值型字段
        numeric_fields = ['total_contacts', 'resolution', 'power_law_slope']
        for field in numeric_fields:
            if not isinstance(data[field], (int, float)):
                error_msg = f"文件 {filename} 字段 {field} 不是数值类型"
                self.validation_errors.append(error_msg)
                self.logger.error(error_msg)
                is_valid = False
        
        # 验证total_contacts为正整数
        if data.get('total_contacts', 0) <= 0:
            warning_msg = f"文件 {filename} total_contacts <= 0"
            self.validation_warnings.append(warning_msg)
            self.logger.warning(warning_msg)
        
        # 验证resolution为正数
        if data.get('resolution', 0) <= 0:
            error_msg = f"文件 {filename} resolution <= 0"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            is_valid = False
        
        return is_valid
    
    def _is_increasing(self, arr: List) -> bool:
        """检查数组是否为递增序列或包含多个递增段（如多染色体数据）"""
        if len(arr) <= 1:
            return True
        
        # 检查是否为单调递增
        is_monotonic = all(arr[i] <= arr[i+1] for i in range(len(arr)-1))
        if is_monotonic:
            return True
        
        # 检查是否为多段递增（如多染色体数据）
        # 每段开始时会有一个大的跳跃（从大值跳到小值）
        segments = []
        current_segment_start = 0
        
        for i in range(len(arr)-1):
            if arr[i] > arr[i+1]:  # 检测到跳跃
                # 验证当前段是否递增
                segment = arr[current_segment_start:i+1]
                if len(segment) > 1 and not all(segment[j] <= segment[j+1] for j in range(len(segment)-1)):
                    return False
                segments.append((current_segment_start, i))
                current_segment_start = i + 1
        
        # 验证最后一段
        if current_segment_start < len(arr):
            segment = arr[current_segment_start:]
            if len(segment) > 1 and not all(segment[j] <= segment[j+1] for j in range(len(segment)-1)):
                return False
            segments.append((current_segment_start, len(arr)-1))
        
        # 如果找到了多个递增段，认为是有效的
        return len(segments) > 0
    
    def _all_non_negative(self, arr: List) -> bool:
        """检查数组是否所有值都非负"""
        return all(x >= 0 for x in arr if x is not None and not np.isnan(x))
    
    def _all_positive(self, arr: List) -> bool:
        """检查数组是否所有值都为正"""
        return all(x > 0 for x in arr if x is not None and not np.isnan(x))
    
    def validate_stage_consistency(self, stage_data: List[Dict], stage: str) -> Tuple[bool, Dict]:
        """
        验证同一stage内数据的一致性
        
        Args:
            stage_data: stage内所有细胞的数据列表
            stage: stage名称
            
        Returns:
            (是否一致, 统计信息)
        """
        if not stage_data:
            error_msg = f"Stage {stage} 没有有效数据"
            self.validation_errors.append(error_msg)
            self.logger.error(error_msg)
            return False, {}
        
        # 收集统计信息
        distances_lengths = []
        resolutions = []
        total_contacts = []
        
        for i, data in enumerate(stage_data):
            distances_lengths.append(len(data.get('distances', [])))
            resolutions.append(data.get('resolution', 0))
            total_contacts.append(data.get('total_contacts', 0))
        
        stats = {
            'cell_count': len(stage_data),
            'distances_length_range': (min(distances_lengths), max(distances_lengths)),
            'distances_length_mode': Counter(distances_lengths).most_common(1)[0] if distances_lengths else (0, 0),
            'resolution_unique': list(set(resolutions)),
            'total_contacts_range': (min(total_contacts), max(total_contacts)),
            'total_contacts_mean': np.mean(total_contacts)
        }
        
        is_consistent = True
        
        # 检查分辨率一致性
        if len(stats['resolution_unique']) > 1:
            warning_msg = f"Stage {stage} 存在不同分辨率: {stats['resolution_unique']}"
            self.validation_warnings.append(warning_msg)
            self.logger.warning(warning_msg)
        
        # 检查distances长度分布
        length_counter = Counter(distances_lengths)
        if len(length_counter) > 3:  # 允许少量长度差异
            warning_msg = f"Stage {stage} distances长度差异较大: {dict(length_counter)}"
            self.validation_warnings.append(warning_msg)
            self.logger.warning(warning_msg)
        
        # 检查总contact数的合理性
        zero_contact_count = sum(1 for x in total_contacts if x == 0)
        if zero_contact_count > 0:
            warning_msg = f"Stage {stage} 有 {zero_contact_count} 个细胞的total_contacts为0"
            self.validation_warnings.append(warning_msg)
            self.logger.warning(warning_msg)
        
        return is_consistent, stats
    
    def align_distances(self, stage_data: List[Dict]) -> Tuple[np.ndarray, np.ndarray]:
        """
        对齐同一stage内所有细胞的distances
        
        Args:
            stage_data: stage内所有细胞的数据列表
            
        Returns:
            (统一的distances数组, 统一的distance_kb数组)
        """
        if not stage_data:
            return np.array([]), np.array([])
        
        # 找到最常见的长度
        distances_lengths = [len(data.get('distances', [])) for data in stage_data]
        length_counter = Counter(distances_lengths)
        target_length = length_counter.most_common(1)[0][0]
        
        # 使用第一个符合target_length的细胞作为模板
        template_data = None
        for data in stage_data:
            if len(data.get('distances', [])) == target_length:
                template_data = data
                break
        
        if template_data is None:
            self.logger.warning(f"无法找到合适的模板数据，使用第一个细胞")
            template_data = stage_data[0]
        
        unified_distances = np.array(template_data['distances'])
        unified_distance_kb = np.array(template_data['distance_kb'])
        
        self.logger.info(f"统一distances长度为: {len(unified_distances)}")
        
        return unified_distances, unified_distance_kb
    
    def create_decay_matrix(self, stage_data: List[Dict], unified_distances: np.ndarray) -> Tuple[np.ndarray, List[str]]:
        """
        创建decay值矩阵
        
        Args:
            stage_data: stage内所有细胞的数据列表
            unified_distances: 统一的distances数组
            
        Returns:
            (细胞×距离的decay值矩阵, 细胞ID列表)
        """
        target_length = len(unified_distances)
        decay_matrix = []
        cell_ids = []
        
        for data in stage_data:
            cell_id = data.get('cell_id', 'unknown')
            decay_values = data.get('decay_values', [])
            
            # 处理长度不匹配
            if len(decay_values) == target_length:
                decay_row = decay_values
            elif len(decay_values) > target_length:
                # 截断
                decay_row = decay_values[:target_length]
                self.logger.warning(f"细胞 {cell_id} 的decay_values被截断从 {len(decay_values)} 到 {target_length}")
            else:
                # 填充NaN
                decay_row = decay_values + [np.nan] * (target_length - len(decay_values))
                self.logger.warning(f"细胞 {cell_id} 的decay_values被填充从 {len(decay_values)} 到 {target_length}")
            
            decay_matrix.append(decay_row)
            cell_ids.append(cell_id)
        
        return np.array(decay_matrix, dtype=np.float32), cell_ids
    
    def get_validation_summary(self) -> Dict:
        """
        获取验证摘要
        
        Returns:
            验证摘要字典
        """
        return {
            'errors': self.validation_errors,
            'warnings': self.validation_warnings,
            'error_count': len(self.validation_errors),
            'warning_count': len(self.validation_warnings),
            'has_errors': len(self.validation_errors) > 0
        }
    
    def log_validation_summary(self) -> None:
        """记录验证摘要到日志"""
        summary = self.get_validation_summary()
        
        if summary['error_count'] > 0:
            self.logger.error(f"验证发现 {summary['error_count']} 个错误:")
            for error in summary['errors']:
                self.logger.error(f"  - {error}")
        
        if summary['warning_count'] > 0:
            self.logger.warning(f"验证发现 {summary['warning_count']} 个警告:")
            for warning in summary['warnings']:
                self.logger.warning(f"  - {warning}")
        
        if not summary['has_errors'] and summary['warning_count'] == 0:
            self.logger.info("数据验证通过，未发现错误或警告")


def test_data_validator():
    """测试数据验证器"""
    # 这里可以添加测试代码
    pass


if __name__ == "__main__":
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    test_data_validator()