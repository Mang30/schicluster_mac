#!/usr/bin/env python3
"""
测试安装和配置脚本
验证所有组件是否正常工作
"""

import sys
import os
from pathlib import Path

def test_imports():
    """测试必需的包导入"""
    print("测试包导入...")
    try:
        import pandas as pd
        print(f"  ✓ pandas {pd.__version__}")
        
        import numpy as np
        print(f"  ✓ numpy {np.__version__}")
        
        import scanpy as sc
        print(f"  ✓ scanpy {sc.__version__}")
        
        import yaml
        print(f"  ✓ pyyaml")
        
        import openpyxl
        print(f"  ✓ openpyxl {openpyxl.__version__}")
        
        import psutil
        print(f"  ✓ psutil {psutil.__version__}")
        
        return True
    except ImportError as e:
        print(f"  ✗ 导入错误: {e}")
        return False

def test_config():
    """测试配置文件"""
    print("\n测试配置文件...")
    try:
        from utils import load_config, ConfigValidator
        
        config = load_config("config.yaml")
        print("  ✓ 配置文件加载成功")
        
        errors = ConfigValidator.validate_config(config)
        if errors:
            print("  ⚠ 配置验证警告:")
            for error in errors:
                print(f"    - {error}")
        else:
            print("  ✓ 配置验证通过")
        
        return len(errors) == 0
    except Exception as e:
        print(f"  ✗ 配置测试失败: {e}")
        return False

def test_paths():
    """测试关键路径"""
    print("\n测试路径...")
    try:
        from utils import load_config
        config = load_config("config.yaml")
        
        paths = config.get('paths', {})
        
        # 测试输入路径
        input_root = paths.get('input_root')
        if input_root and Path(input_root).exists():
            print(f"  ✓ 输入目录存在: {input_root}")
        else:
            print(f"  ⚠ 输入目录不存在: {input_root}")
        
        # 测试输出路径
        output_dir = paths.get('output_dir')
        if output_dir:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            print(f"  ✓ 输出目录已创建: {output_dir}")
        
        # 测试日志路径
        log_dir = paths.get('log_dir')
        if log_dir:
            Path(log_dir).mkdir(parents=True, exist_ok=True)
            print(f"  ✓ 日志目录已创建: {log_dir}")
        
        # 测试元数据文件
        metadata_file = paths.get('metadata_file')
        if metadata_file and Path(metadata_file).exists():
            print(f"  ✓ 元数据文件存在: {metadata_file}")
        else:
            print(f"  ⚠ 元数据文件不存在: {metadata_file}")
            print("    请确保 GSE223917_HiRES_emb_metadata.xlsx 文件在正确位置")
        
        return True
    except Exception as e:
        print(f"  ✗ 路径测试失败: {e}")
        return False

def test_modules():
    """测试本地模块"""
    print("\n测试本地模块...")
    try:
        from metadata_processor import MetadataProcessor
        print("  ✓ metadata_processor 模块正常")
        
        from data_validator import DataValidator
        print("  ✓ data_validator 模块正常")
        
        from utils import setup_logging, ProgressTracker, MemoryMonitor
        print("  ✓ utils 模块正常")
        
        return True
    except ImportError as e:
        print(f"  ✗ 模块导入错误: {e}")
        return False

def test_sample_data():
    """测试样本数据"""
    print("\n测试样本数据...")
    try:
        from utils import load_config, get_file_list
        config = load_config("config.yaml")
        
        input_root = config['paths']['input_root']
        
        # 检查是否有 stage 目录
        for stage in config['stages']:
            stage_dir = Path(input_root) / f"stage_{stage}"
            if stage_dir.exists():
                json_files = get_file_list(str(stage_dir), "*_decay_profile.json")
                if json_files:
                    print(f"  ✓ Stage {stage}: 找到 {len(json_files)} 个 JSON 文件")
                else:
                    print(f"  ⚠ Stage {stage}: 未找到 JSON 文件")
            else:
                print(f"  ⚠ Stage {stage}: 目录不存在 ({stage_dir})")
        
        return True
    except Exception as e:
        print(f"  ✗ 样本数据测试失败: {e}")
        return False

def test_write_permissions():
    """测试写入权限"""
    print("\n测试写入权限...")
    try:
        from utils import load_config
        config = load_config("config.yaml")
        
        # 测试输出目录写入权限
        output_dir = Path(config['paths']['output_dir'])
        test_file = output_dir / "test_write.tmp"
        
        with open(test_file, 'w') as f:
            f.write("test")
        
        test_file.unlink()  # 删除测试文件
        print("  ✓ 输出目录可写")
        
        # 测试日志目录写入权限
        log_dir = Path(config['paths']['log_dir'])
        test_log = log_dir / "test_log.tmp"
        
        with open(test_log, 'w') as f:
            f.write("test log")
        
        test_log.unlink()  # 删除测试文件
        print("  ✓ 日志目录可写")
        
        return True
    except Exception as e:
        print(f"  ✗ 写入权限测试失败: {e}")
        return False

def main():
    """主测试函数"""
    print("JSON to H5AD 转换工具 - 安装测试")
    print("=" * 50)
    
    tests = [
        ("包导入", test_imports),
        ("配置文件", test_config),
        ("路径检查", test_paths),
        ("本地模块", test_modules),
        ("样本数据", test_sample_data),
        ("写入权限", test_write_permissions)
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"  ✗ {test_name} 测试异常: {e}")
            results.append((test_name, False))
    
    print("\n" + "=" * 50)
    print("测试结果汇总:")
    
    all_passed = True
    for test_name, result in results:
        status = "✓ 通过" if result else "✗ 失败"
        print(f"  {test_name}: {status}")
        if not result:
            all_passed = False
    
    print("\n" + "=" * 50)
    if all_passed:
        print("🎉 所有测试通过! 系统已准备就绪。")
        print("\n下一步:")
        print("1. 确保元数据文件 GSE223917_HiRES_emb_metadata.xlsx 在正确位置")
        print("2. 运行单个 stage 测试:")
        print("   python convert_json_to_h5ad.py --stage E75 --debug")
        print("3. 或运行批处理:")
        print("   ./process_all_stages.sh --debug")
    else:
        print("⚠ 部分测试失败，请检查上述错误并修复。")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())