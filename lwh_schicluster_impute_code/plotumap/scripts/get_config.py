#!/usr/bin/env python3
"""
路径辅助脚本 - 为Shell脚本提供路径配置支持
"""

import sys
import os

# 添加core目录到Python路径
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
core_dir = os.path.join(project_root, 'core')
sys.path.insert(0, core_dir)

try:
    from config import CONFIG
    
    if len(sys.argv) > 1:
        attr_name = sys.argv[1]
        if hasattr(CONFIG, attr_name):
            print(getattr(CONFIG, attr_name))
        else:
            print(f"ERROR: Unknown config attribute: {attr_name}", file=sys.stderr)
            sys.exit(1)
    else:
        # 打印所有可用的配置
        print("Available configurations:")
        for attr in dir(CONFIG):
            if not attr.startswith('_') and attr.isupper():
                print(f"  {attr}: {getattr(CONFIG, attr)}")

except ImportError as e:
    print(f"ERROR: Failed to import config: {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
