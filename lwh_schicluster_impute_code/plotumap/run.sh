#!/bin/bash

# 🚀 快速启动脚本 - 调用主处理脚本
# 使用相对路径，自动定位到 scripts/run.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_SCRIPT="$SCRIPT_DIR/scripts/run.sh"

# 检查主脚本是否存在
if [[ ! -f "$MAIN_SCRIPT" ]]; then
    echo "❌ 错误: 找不到主脚本 $MAIN_SCRIPT"
    exit 1
fi

# 调用主脚本，传递所有参数
exec "$MAIN_SCRIPT" "$@"
