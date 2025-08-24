#!/usr/bin/env python3
"""
配置文件 - 管理所有文件路径
"""

import os
from pathlib import Path

# 获取项目根目录
PROJECT_ROOT = Path(__file__).parent.parent  # plotumap 目录
BASE_DIR = PROJECT_ROOT.parent.parent  # 0_HiRES 目录

class PathConfig:
    """路径配置类"""
    
    # 基础目录
    BASE = str(BASE_DIR)
    PROJECT_ROOT_STR = str(PROJECT_ROOT)
    
    # 数据目录
    DATA_DIR = str(BASE_DIR / "data")
    OUTPUT_DIR = str(BASE_DIR / "output")
    MERGED_MATRICES = str(BASE_DIR / "output" / "merged_matrices")
    IMPUTED_MATRICES = str(BASE_DIR / "output" / "imputed_matrices_by_stage")
    H5AD_OUTPUT = str(BASE_DIR / "output" / "hic_h5ad_files")
    
    # 配置文件
    CONFIG_DIR = str(PROJECT_ROOT / "config")
    STAGE_FILES_MAPPING = str(PROJECT_ROOT / "config" / "stage_files_mapping.csv")
    COLOR_MAPPING = str(PROJECT_ROOT / "config" / "color_mapping.json")
    
    # 参考文件
    CHROM_SIZES = str(BASE_DIR / "mm10_chrom_sizes.txt")
    
    # 文档和资源
    DOCS_DIR = str(PROJECT_ROOT / "docs")
    
    @classmethod
    def get_stage_dir(cls, stage: str) -> str:
        """获取指定阶段的目录路径"""
        return str(Path(cls.IMPUTED_MATRICES) / stage)
    
    @classmethod
    def get_output_path(cls, stage: str, suffix: str = "") -> str:
        """获取输出文件路径"""
        filename = f"{stage}_hic{suffix}.h5ad"
        return str(Path(cls.H5AD_OUTPUT) / filename)
    
    @classmethod
    def check_paths(cls) -> dict:
        """检查所有路径是否存在"""
        paths_status = {}
        
        # 检查目录
        for attr_name in dir(cls):
            if attr_name.isupper() and attr_name.endswith('_DIR'):
                path = getattr(cls, attr_name)
                paths_status[attr_name] = os.path.exists(path)
        
        # 检查文件
        file_paths = {
            'STAGE_FILES_MAPPING': cls.STAGE_FILES_MAPPING,
            'COLOR_MAPPING': cls.COLOR_MAPPING,
            'CHROM_SIZES': cls.CHROM_SIZES
        }
        
        for name, path in file_paths.items():
            paths_status[name] = os.path.exists(path)
            
        return paths_status

# 创建全局配置实例
CONFIG = PathConfig()

if __name__ == '__main__':
    print("📁 路径配置检查:")
    print(f"项目根目录: {CONFIG.PROJECT_ROOT_STR}")
    print(f"基础目录: {CONFIG.BASE}")
    print("")
    
    status = CONFIG.check_paths()
    for name, exists in status.items():
        status_icon = "✅" if exists else "❌"
        path = getattr(CONFIG, name)
        print(f"{status_icon} {name}: {path}")
