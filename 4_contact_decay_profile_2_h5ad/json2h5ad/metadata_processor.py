#!/usr/bin/env python3
"""
元数据处理模块
用于读取和处理Excel元数据文件，实现细胞ID匹配和字段映射
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, List
import warnings

class MetadataProcessor:
    """元数据处理器"""
    
    def __init__(self, metadata_file: str):
        """
        初始化元数据处理器
        
        Args:
            metadata_file: Excel元数据文件路径
        """
        self.metadata_file = Path(metadata_file)
        self.metadata_df = None
        self.logger = logging.getLogger(__name__)
        
        # 字段映射配置
        self.field_mapping = {
            'Cellname': 'cell_id',
            'Stage': 'stage',
            'Celltype': 'celltype', 
            'Cellcycle phase': 'cellcycle_phase',
            'CDPS cluster': 'cdps_cluster',
            'Sub_k_cluster': 'sub_k_cluster',
            'G1S.Score': 'g1s_score',
            'G2M.Score': 'g2m_score',
            'repli score': 'repli_score',
            'Raw contacts': 'raw_contacts',
            'Clean3 contacts': 'clean_contacts'
        }
        
        self._load_metadata()
    
    def _load_metadata(self) -> None:
        """加载Excel元数据文件"""
        try:
            if not self.metadata_file.exists():
                raise FileNotFoundError(f"元数据文件不存在: {self.metadata_file}")
            
            self.logger.info(f"加载元数据文件: {self.metadata_file}")
            self.metadata_df = pd.read_excel(self.metadata_file)
            
            # 检查必需字段
            required_fields = ['Cellname', 'Stage', 'Celltype']
            missing_fields = [field for field in required_fields 
                            if field not in self.metadata_df.columns]
            if missing_fields:
                raise ValueError(f"元数据文件缺少必需字段: {missing_fields}")
            
            # 清理数据
            self._clean_metadata()
            
            self.logger.info(f"成功加载元数据: {len(self.metadata_df)} 个细胞记录")
            
        except Exception as e:
            self.logger.error(f"加载元数据文件失败: {e}")
            raise
    
    def _clean_metadata(self) -> None:
        """清理和标准化元数据"""
        # 去除空白字符
        string_cols = self.metadata_df.select_dtypes(include=['object']).columns
        for col in string_cols:
            if col in self.metadata_df.columns:
                self.metadata_df[col] = self.metadata_df[col].astype(str).str.strip()
        
        # 标准化细胞名称 (去除可能的前后空白)
        self.metadata_df['Cellname'] = self.metadata_df['Cellname'].str.strip()
        
        # 标准化stage格式
        if 'Stage' in self.metadata_df.columns:
            self.metadata_df['Stage'] = self.metadata_df['Stage'].str.upper()
        
        # 处理数值型字段
        numeric_fields = ['G1S.Score', 'G2M.Score', 'repli score', 'Raw contacts', 'Clean3 contacts']
        for field in numeric_fields:
            if field in self.metadata_df.columns:
                self.metadata_df[field] = pd.to_numeric(self.metadata_df[field], errors='coerce')
    
    def extract_cell_id_from_filename(self, filename: str) -> str:
        """
        从JSON文件名中提取细胞ID
        
        Args:
            filename: JSON文件名 (如 'GasdE701023_decay_profile.json')
            
        Returns:
            细胞ID (如 'GasdE701023')
        """
        # 移除路径和扩展名
        basename = Path(filename).stem
        
        # 移除 '_decay_profile' 后缀
        if basename.endswith('_decay_profile'):
            cell_id = basename.replace('_decay_profile', '')
        else:
            cell_id = basename
        
        return cell_id
    
    def get_cell_metadata(self, cell_id: str) -> Dict:
        """
        根据细胞ID获取元数据
        
        Args:
            cell_id: 细胞ID
            
        Returns:
            包含元数据的字典
        """
        # 查找匹配的记录
        mask = self.metadata_df['Cellname'] == cell_id
        matched_rows = self.metadata_df[mask]
        
        if matched_rows.empty:
            self.logger.warning(f"未找到细胞ID '{cell_id}' 的元数据，使用默认值")
            return self._get_default_metadata(cell_id)
        
        if len(matched_rows) > 1:
            self.logger.warning(f"细胞ID '{cell_id}' 有多条记录，使用第一条")
        
        # 获取第一条匹配记录
        row = matched_rows.iloc[0]
        
        # 构建元数据字典
        metadata = {}
        for excel_field, h5ad_field in self.field_mapping.items():
            if excel_field in row.index:
                value = row[excel_field]
                # 处理NaN值
                if pd.isna(value):
                    if excel_field in ['G1S.Score', 'G2M.Score', 'repli score', 'Raw contacts', 'Clean3 contacts']:
                        metadata[h5ad_field] = np.nan
                    else:
                        metadata[h5ad_field] = 'Unknown'
                else:
                    metadata[h5ad_field] = value
            else:
                metadata[h5ad_field] = 'Unknown' if h5ad_field not in ['g1s_score', 'g2m_score', 'repli_score', 'raw_contacts', 'clean_contacts'] else np.nan
        
        return metadata
    
    def _get_default_metadata(self, cell_id: str) -> Dict:
        """
        获取缺失细胞的默认元数据
        
        Args:
            cell_id: 细胞ID
            
        Returns:
            默认元数据字典
        """
        return {
            'cell_id': cell_id,
            'stage': 'Unknown',
            'celltype': 'Unknown',
            'cellcycle_phase': 'Unknown',
            'cdps_cluster': 'Unknown',
            'sub_k_cluster': 'Unknown',
            'g1s_score': np.nan,
            'g2m_score': np.nan,
            'repli_score': np.nan,
            'raw_contacts': np.nan,
            'clean_contacts': np.nan
        }
    
    def get_stage_cells(self, stage: str) -> List[str]:
        """
        获取指定stage的所有细胞ID
        
        Args:
            stage: 发育阶段 (如 'E70')
            
        Returns:
            细胞ID列表
        """
        if self.metadata_df is None:
            return []
        
        stage_upper = stage.upper()
        mask = self.metadata_df['Stage'] == stage_upper
        stage_cells = self.metadata_df[mask]['Cellname'].tolist()
        
        self.logger.info(f"Stage {stage} 包含 {len(stage_cells)} 个细胞")
        return stage_cells
    
    def get_celltype_counts(self, stage: Optional[str] = None) -> Dict[str, int]:
        """
        获取细胞类型分布统计
        
        Args:
            stage: 可选的stage过滤
            
        Returns:
            细胞类型计数字典
        """
        if self.metadata_df is None:
            return {}
        
        df = self.metadata_df
        if stage:
            stage_upper = stage.upper()
            df = df[df['Stage'] == stage_upper]
        
        if 'Celltype' in df.columns:
            counts = df['Celltype'].value_counts().to_dict()
        else:
            counts = {}
        
        return counts
    
    def get_stage_summary(self) -> Dict[str, Dict]:
        """
        获取各stage的细胞数量统计
        
        Returns:
            stage统计信息字典
        """
        if self.metadata_df is None:
            return {}
        
        summary = {}
        if 'Stage' in self.metadata_df.columns:
            for stage in self.metadata_df['Stage'].unique():
                if pd.notna(stage):
                    stage_df = self.metadata_df[self.metadata_df['Stage'] == stage]
                    summary[stage] = {
                        'cell_count': len(stage_df),
                        'celltype_counts': self.get_celltype_counts(stage)
                    }
        
        return summary
    
    def validate_metadata(self) -> Tuple[bool, List[str]]:
        """
        验证元数据完整性
        
        Returns:
            (是否有效, 错误信息列表)
        """
        errors = []
        
        if self.metadata_df is None:
            errors.append("元数据未加载")
            return False, errors
        
        # 检查必需字段
        required_fields = ['Cellname', 'Stage', 'Celltype']
        for field in required_fields:
            if field not in self.metadata_df.columns:
                errors.append(f"缺少必需字段: {field}")
            elif self.metadata_df[field].isna().any():
                na_count = self.metadata_df[field].isna().sum()
                errors.append(f"字段 {field} 有 {na_count} 个空值")
        
        # 检查细胞ID唯一性
        if 'Cellname' in self.metadata_df.columns:
            duplicated = self.metadata_df['Cellname'].duplicated().any()
            if duplicated:
                dup_count = self.metadata_df['Cellname'].duplicated().sum()
                errors.append(f"发现 {dup_count} 个重复的细胞ID")
        
        # 检查stage格式
        if 'Stage' in self.metadata_df.columns:
            valid_stages = {'E70', 'E75', 'E80', 'E85', 'E95', 'EX05', 'EX15'}
            invalid_stages = set(self.metadata_df['Stage'].unique()) - valid_stages - {np.nan}
            if invalid_stages:
                errors.append(f"发现无效的stage: {invalid_stages}")
        
        is_valid = len(errors) == 0
        return is_valid, errors


def test_metadata_processor():
    """测试元数据处理器"""
    # 这里可以添加测试代码
    pass


if __name__ == "__main__":
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    test_metadata_processor()