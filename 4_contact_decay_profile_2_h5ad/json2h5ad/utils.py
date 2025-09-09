#!/usr/bin/env python3
"""
工具函数和日志记录系统
提供通用工具函数、日志配置和性能监控
"""

import os
import sys
import logging
import time
import psutil
import yaml
import json
from pathlib import Path
from typing import Dict, Any, Optional, List
from datetime import datetime
import numpy as np
import pandas as pd
from contextlib import contextmanager


def load_config(config_file: str = "config.yaml") -> Dict[str, Any]:
    """
    加载配置文件
    
    Args:
        config_file: 配置文件路径
        
    Returns:
        配置字典
    """
    config_path = Path(config_file)
    if not config_path.exists():
        raise FileNotFoundError(f"配置文件不存在: {config_file}")
    
    with open(config_file, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    return config


def setup_logging(config: Dict[str, Any], stage: Optional[str] = None) -> logging.Logger:
    """
    配置日志系统
    
    Args:
        config: 配置字典
        stage: 可选的stage名称，用于创建特定的日志文件
        
    Returns:
        配置好的logger
    """
    log_config = config.get('logging', {})
    
    # 创建日志目录
    log_dir = Path(config['paths']['log_dir'])
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # 设置日志级别
    log_level = getattr(logging, log_config.get('level', 'INFO').upper())
    
    # 配置根logger
    logger = logging.getLogger()
    logger.setLevel(log_level)
    
    # 清除已有的handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 创建formatter
    formatter = logging.Formatter(
        log_config.get('format', '%(asctime)s - %(name)s - %(levelname)s - %(message)s'),
        datefmt=log_config.get('date_format', '%Y-%m-%d %H:%M:%S')
    )
    
    # 控制台handler
    if log_config.get('console_output', True):
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    # 文件handler
    if stage and log_config.get('stage_specific_logs', True):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_filename = log_config.get('log_filename_format', 'stage_{stage}_processing_{timestamp}.log').format(
            stage=stage, timestamp=timestamp
        )
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_filename = f"json2h5ad_processing_{timestamp}.log"
    
    log_file = log_dir / log_filename
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    logger.info(f"日志系统初始化完成，日志文件: {log_file}")
    return logger


class ProgressTracker:
    """进度跟踪器"""
    
    def __init__(self, total: int, description: str = "Processing", update_interval: int = 10):
        """
        初始化进度跟踪器
        
        Args:
            total: 总任务数
            description: 任务描述
            update_interval: 更新间隔
        """
        self.total = total
        self.description = description
        self.update_interval = update_interval
        self.current = 0
        self.start_time = time.time()
        self.last_update = 0
        self.logger = logging.getLogger(__name__)
    
    def update(self, count: int = 1) -> None:
        """
        更新进度
        
        Args:
            count: 增加的任务数
        """
        self.current += count
        
        if (self.current - self.last_update >= self.update_interval or 
            self.current == self.total):
            
            progress = self.current / self.total * 100
            elapsed = time.time() - self.start_time
            
            if self.current > 0:
                eta = elapsed * (self.total - self.current) / self.current
                eta_str = format_duration(eta)
            else:
                eta_str = "unknown"
            
            self.logger.info(f"{self.description}: {self.current}/{self.total} "
                           f"({progress:.1f}%) - 已用时: {format_duration(elapsed)} - "
                           f"预计剩余: {eta_str}")
            
            self.last_update = self.current
    
    def finish(self) -> None:
        """完成进度跟踪"""
        total_time = time.time() - self.start_time
        self.logger.info(f"{self.description} 完成! 总用时: {format_duration(total_time)}")


class MemoryMonitor:
    """内存监控器"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        初始化内存监控器
        
        Args:
            logger: 日志记录器
        """
        self.logger = logger or logging.getLogger(__name__)
        self.process = psutil.Process(os.getpid())
        self.start_memory = self.get_memory_usage()
    
    def get_memory_usage(self) -> float:
        """
        获取当前内存使用量 (MB)
        
        Returns:
            内存使用量 (MB)
        """
        return self.process.memory_info().rss / 1024 / 1024
    
    def log_memory_usage(self, description: str = "Memory usage") -> None:
        """
        记录内存使用情况
        
        Args:
            description: 描述信息
        """
        current_memory = self.get_memory_usage()
        memory_increase = current_memory - self.start_memory
        
        self.logger.info(f"{description}: {current_memory:.1f} MB "
                        f"(+{memory_increase:.1f} MB from start)")
    
    @contextmanager
    def monitor_block(self, description: str = "Code block"):
        """
        监控代码块的内存使用
        
        Args:
            description: 代码块描述
        """
        start_memory = self.get_memory_usage()
        start_time = time.time()
        
        try:
            yield
        finally:
            end_memory = self.get_memory_usage()
            end_time = time.time()
            
            memory_change = end_memory - start_memory
            duration = end_time - start_time
            
            self.logger.info(f"{description} 完成 - 用时: {format_duration(duration)} - "
                           f"内存变化: {memory_change:+.1f} MB")


def format_duration(seconds: float) -> str:
    """
    格式化持续时间
    
    Args:
        seconds: 秒数
        
    Returns:
        格式化的时间字符串
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


def format_file_size(size_bytes: int) -> str:
    """
    格式化文件大小
    
    Args:
        size_bytes: 字节数
        
    Returns:
        格式化的文件大小字符串
    """
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


def create_directory(directory: str, exist_ok: bool = True) -> Path:
    """
    创建目录
    
    Args:
        directory: 目录路径
        exist_ok: 如果目录已存在是否报错
        
    Returns:
        Path对象
    """
    dir_path = Path(directory)
    dir_path.mkdir(parents=True, exist_ok=exist_ok)
    return dir_path


def backup_file(file_path: str, backup_suffix: str = ".backup") -> Optional[Path]:
    """
    备份文件
    
    Args:
        file_path: 源文件路径
        backup_suffix: 备份文件后缀
        
    Returns:
        备份文件路径，如果源文件不存在则返回None
    """
    source_path = Path(file_path)
    if not source_path.exists():
        return None
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = source_path.with_suffix(f"{backup_suffix}_{timestamp}")
    
    import shutil
    shutil.copy2(source_path, backup_path)
    
    return backup_path


def validate_environment(config: Dict[str, Any]) -> List[str]:
    """
    验证运行环境
    
    Args:
        config: 配置字典
        
    Returns:
        错误信息列表
    """
    errors = []
    
    # 检查Python路径
    python_path = config.get('environment', {}).get('python_path')
    if python_path:
        if not Path(python_path).exists():
            errors.append(f"Python路径不存在: {python_path}")
    
    # 检查必需的包
    required_packages = config.get('environment', {}).get('required_packages', {})
    for package, min_version in required_packages.items():
        try:
            __import__(package)
        except ImportError:
            errors.append(f"缺少必需的包: {package}")
    
    # 检查路径
    paths = config.get('paths', {})
    for path_name, path_value in paths.items():
        if path_name in ['input_root', 'metadata_file']:  # 检查输入路径
            if not Path(path_value).exists():
                errors.append(f"路径不存在: {path_name} = {path_value}")
    
    return errors


def save_json(data: Dict[str, Any], file_path: str, indent: int = 2) -> None:
    """
    保存JSON文件
    
    Args:
        data: 要保存的数据
        file_path: 文件路径
        indent: 缩进
    """
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=indent, default=str)


def load_json(file_path: str) -> Dict[str, Any]:
    """
    加载JSON文件
    
    Args:
        file_path: 文件路径
        
    Returns:
        数据字典
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def check_disk_space(path: str, required_gb: float = 1.0) -> bool:
    """
    检查磁盘空间
    
    Args:
        path: 检查路径
        required_gb: 需要的空间 (GB)
        
    Returns:
        是否有足够空间
    """
    try:
        statvfs = os.statvfs(path)
        available_bytes = statvfs.f_frsize * statvfs.f_bavail
        available_gb = available_bytes / (1024 ** 3)
        return available_gb >= required_gb
    except:
        return True  # 如果检查失败，假设有足够空间


def get_file_list(directory: str, pattern: str = "*_decay_profile.json") -> List[Path]:
    """
    获取目录中匹配模式的文件列表
    
    Args:
        directory: 目录路径
        pattern: 文件模式
        
    Returns:
        文件路径列表
    """
    dir_path = Path(directory)
    if not dir_path.exists():
        return []
    
    return sorted(dir_path.glob(pattern))


class ConfigValidator:
    """配置验证器"""
    
    @staticmethod
    def validate_config(config: Dict[str, Any]) -> List[str]:
        """
        验证配置文件
        
        Args:
            config: 配置字典
            
        Returns:
            错误信息列表
        """
        errors = []
        
        # 检查必需的顶级字段
        required_sections = ['paths', 'stages', 'processing']
        for section in required_sections:
            if section not in config:
                errors.append(f"配置文件缺少必需部分: {section}")
        
        # 检查路径配置
        if 'paths' in config:
            required_paths = ['input_root', 'output_dir', 'metadata_file', 'log_dir']
            for path_name in required_paths:
                if path_name not in config['paths']:
                    errors.append(f"路径配置缺少: {path_name}")
        
        # 检查stages配置
        if 'stages' in config:
            valid_stages = {'E70', 'E75', 'E80', 'E85', 'E95', 'EX05', 'EX15'}
            configured_stages = set(config['stages'])
            invalid_stages = configured_stages - valid_stages
            if invalid_stages:
                errors.append(f"配置了无效的stages: {invalid_stages}")
        
        return errors


def get_system_info() -> Dict[str, Any]:
    """
    获取系统信息
    
    Returns:
        系统信息字典
    """
    return {
        'platform': sys.platform,
        'python_version': sys.version,
        'cpu_count': psutil.cpu_count(),
        'memory_total_gb': psutil.virtual_memory().total / (1024 ** 3),
        'memory_available_gb': psutil.virtual_memory().available / (1024 ** 3),
        'current_time': datetime.now().isoformat()
    }


if __name__ == "__main__":
    # 测试代码
    config = load_config()
    logger = setup_logging(config)
    logger.info("工具模块测试完成")