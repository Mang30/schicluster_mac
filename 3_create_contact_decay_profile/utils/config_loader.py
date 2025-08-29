#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
配置文件加载工具

功能：
1. 加载YAML配置文件
2. 提供配置验证
3. 支持配置合并和覆盖

作者：Claude Code Assistant
日期：2025-08-29
"""

import os
import sys
import yaml
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List

class ConfigLoader:
    """配置加载器"""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        初始化配置加载器
        
        Args:
            config_path: 配置文件路径，如果为None则使用默认路径
        """
        if config_path is None:
            # 使用默认配置路径
            script_dir = Path(__file__).parent
            config_path = script_dir.parent / 'configs' / 'analysis_config.yaml'
        
        self.config_path = Path(config_path)
        self.config = {}
        self.logger = logging.getLogger(__name__)
        
    def load_config(self) -> Dict[str, Any]:
        """
        加载配置文件
        
        Returns:
            Dict: 配置字典
        """
        try:
            if not self.config_path.exists():
                raise FileNotFoundError(f"配置文件不存在: {self.config_path}")
            
            with open(self.config_path, 'r', encoding='utf-8') as f:
                self.config = yaml.safe_load(f)
            
            self.logger.info(f"成功加载配置文件: {self.config_path}")
            return self.config
            
        except Exception as e:
            self.logger.error(f"加载配置文件失败: {e}")
            return {}
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        获取配置值，支持点号分隔的嵌套键
        
        Args:
            key: 配置键，支持 'section.subsection.key' 格式
            default: 默认值
            
        Returns:
            配置值
        """
        if not self.config:
            self.load_config()
        
        try:
            keys = key.split('.')
            value = self.config
            
            for k in keys:
                if isinstance(value, dict) and k in value:
                    value = value[k]
                else:
                    return default
            
            return value
            
        except Exception:
            return default
    
    def get_stage_config(self) -> Dict[str, Any]:
        """获取发育阶段相关配置"""
        return self.get('stages', {})
    
    def get_analysis_config(self) -> Dict[str, Any]:
        """获取分析参数配置"""
        return self.get('analysis', {})
    
    def get_parallel_config(self) -> Dict[str, Any]:
        """获取并行处理配置"""
        return self.get('parallel', {})
    
    def get_output_config(self) -> Dict[str, Any]:
        """获取输出配置"""
        return self.get('output', {})
    
    def get_paths_config(self) -> Dict[str, Any]:
        """获取路径配置"""
        return self.get('paths', {})
    
    def validate_config(self) -> List[str]:
        """
        验证配置文件
        
        Returns:
            List[str]: 错误信息列表，空列表表示验证通过
        """
        errors = []
        
        if not self.config:
            errors.append("配置为空")
            return errors
        
        # 检查必需的配置节
        required_sections = ['paths', 'stages', 'analysis', 'output']
        for section in required_sections:
            if section not in self.config:
                errors.append(f"缺少必需的配置节: {section}")
        
        # 检查路径配置
        paths = self.get('paths', {})
        if 'data_root' not in paths:
            errors.append("缺少数据根目录配置: paths.data_root")
        if 'output_root' not in paths:
            errors.append("缺少输出根目录配置: paths.output_root")
        
        # 检查发育阶段配置
        stages = self.get('stages', {})
        if 'stage_list' not in stages or not stages['stage_list']:
            errors.append("缺少发育阶段列表配置: stages.stage_list")
        
        # 检查分析配置
        analysis = self.get('analysis', {})
        if 'decay_profile' not in analysis:
            errors.append("缺少衰减曲线分析配置: analysis.decay_profile")
        
        return errors
    
    def create_directories(self) -> bool:
        """
        根据配置创建必要的目录
        
        Returns:
            bool: 是否成功创建所有目录
        """
        try:
            paths = self.get_paths_config()
            
            # 需要创建的目录列表
            dirs_to_create = [
                'output_root',
                'log_dir', 
                'temp_dir',
                'plots_dir'
            ]
            
            for dir_key in dirs_to_create:
                dir_path = paths.get(dir_key)
                if dir_path:
                    Path(dir_path).mkdir(parents=True, exist_ok=True)
                    self.logger.info(f"创建目录: {dir_path}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"创建目录失败: {e}")
            return False
    
    def get_stage_colors(self) -> Dict[str, str]:
        """获取发育阶段颜色映射"""
        return self.get('stages.stage_colors', {})
    
    def get_stage_descriptions(self) -> Dict[str, str]:
        """获取发育阶段描述映射"""
        return self.get('stages.stage_descriptions', {})
    
    def override_config(self, overrides: Dict[str, Any]) -> None:
        """
        覆盖配置值
        
        Args:
            overrides: 要覆盖的配置字典，支持点号分隔的键
        """
        for key, value in overrides.items():
            self._set_nested_value(key, value)
    
    def _set_nested_value(self, key: str, value: Any) -> None:
        """设置嵌套配置值"""
        if not self.config:
            self.load_config()
        
        keys = key.split('.')
        config_ref = self.config
        
        # 导航到目标位置
        for k in keys[:-1]:
            if k not in config_ref:
                config_ref[k] = {}
            config_ref = config_ref[k]
        
        # 设置值
        config_ref[keys[-1]] = value
    
    def save_config(self, output_path: Optional[str] = None) -> bool:
        """
        保存配置到文件
        
        Args:
            output_path: 输出路径，如果为None则覆盖原文件
            
        Returns:
            bool: 是否保存成功
        """
        try:
            save_path = Path(output_path) if output_path else self.config_path
            
            with open(save_path, 'w', encoding='utf-8') as f:
                yaml.dump(self.config, f, default_flow_style=False, 
                         allow_unicode=True, indent=2)
            
            self.logger.info(f"配置已保存到: {save_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"保存配置失败: {e}")
            return False


# 全局配置实例
_config_instance = None

def get_config(config_path: Optional[str] = None) -> ConfigLoader:
    """
    获取全局配置实例
    
    Args:
        config_path: 配置文件路径
        
    Returns:
        ConfigLoader: 配置加载器实例
    """
    global _config_instance
    
    if _config_instance is None:
        _config_instance = ConfigLoader(config_path)
        _config_instance.load_config()
    
    return _config_instance


def main():
    """测试函数"""
    import logging
    logging.basicConfig(level=logging.INFO)
    
    # 创建配置加载器
    config_loader = ConfigLoader()
    
    # 加载配置
    config = config_loader.load_config()
    
    if config:
        print("✅ 配置加载成功")
        
        # 验证配置
        errors = config_loader.validate_config()
        if errors:
            print("❌ 配置验证失败:")
            for error in errors:
                print(f"  - {error}")
        else:
            print("✅ 配置验证通过")
        
        # 显示一些关键配置
        print(f"\n📊 项目名称: {config_loader.get('project_name', 'Unknown')}")
        print(f"🔬 发育阶段: {config_loader.get('stages.stage_list', [])}")
        print(f"⚡ 并行处理: {config_loader.get('parallel.enable_parallel', False)}")
        print(f"📈 最大分析距离: {config_loader.get('analysis.decay_profile.max_distance', 'Unknown')}")
        
        # 创建目录
        if config_loader.create_directories():
            print("✅ 目录创建成功")
        else:
            print("❌ 目录创建失败")
    
    else:
        print("❌ 配置加载失败")


if __name__ == "__main__":
    main()