# E75阶段全量处理脚本使用说明

## 脚本位置
脚本位于: `/home/duxuyan/Projects/schicluster_mac/3_create_contact_decay_profile/scripts/run_stage_E75.sh`

## 推荐执行方式

### 方法1: 从脚本所在目录执行（推荐）
```bash
cd /home/duxuyan/Projects/schicluster_mac/3_create_contact_decay_profile/scripts
./run_stage_E75.sh
```

### 方法2: 使用绝对路径执行
```bash
/home/duxuyan/Projects/schicluster_mac/3_create_contact_decay_profile/scripts/run_stage_E75.sh
```

### 方法3: 从项目根目录执行
```bash
cd /home/duxuyan/Projects/schicluster_mac/3_create_contact_decay_profile
scripts/run_stage_E75.sh
```

## 命令行参数

- `-w, --max-workers N`: 设置并行处理的worker数量（默认: 4）
- `-d, --max-distance N`: 设置最大分析距离（默认: 1000）
- `-n, --max-files N`: 限制处理的最大文件数（默认: 不限制，全量处理）
- `--data-root PATH`: 指定数据根目录路径
- `-h, --help`: 显示帮助信息

## 示例

```bash
# 基本执行（全量处理）
./run_stage_E75.sh

# 设置并行worker数为8，最大分析距离为2000
./run_stage_E75.sh -w 8 -d 2000

# 限制只处理前10个文件
./run_stage_E75.sh -n 10

# 指定数据目录路径
./run_stage_E75.sh --data-root /path/to/your/data
```

## 输出目录
处理结果将保存在: `/home/duxuyan/Projects/schicluster_mac/3_create_contact_decay_profile/outputs/stage_E75/`

## 注意事项

1. 确保已安装micromamba并配置好`3_schicluster_python38`环境
2. 脚本具有路径自适应能力，可以从不同目录执行
3. 如果数据目录无法自动找到，脚本会提示手动输入路径
4. 脚本会自动检查所有依赖项并在需要时安装