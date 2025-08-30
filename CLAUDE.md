# HiRES 项目 - Claude 使用指南

## 环境配置

### macOS 下使用 micromamba 环境

本项目使用 `schicluster` 环境进行数据处理和分析。

#### 正确的 Python 命令
```bash
# 直接使用环境中的 Python（推荐）
/Users/wuhaoliu/mamba/envs/schicluster/bin/python script.py

# 或者使用 micromamba run（如果配置正确）
/opt/homebrew/bin/micromamba run -n schicluster python script.py
```

#### 环境路径
- **micromamba 可执行文件**: `/opt/homebrew/bin/micromamba`
- **schicluster 环境路径**: `/Users/wuhaoliu/mamba/envs/schicluster`
- **环境中的 Python**: `/Users/wuhaoliu/mamba/envs/schicluster/bin/python`

#### 检查环境
```bash
# 列出所有环境
/opt/homebrew/bin/micromamba env list

# 检查 pandas 是否可用
/Users/wuhaoliu/mamba/envs/schicluster/bin/python -c "import pandas; print('pandas version:', pandas.__version__)"
```

## 项目数据路径

### Hi-C 数据文件位置
- **数据目录**: `/Volumes/SumSung500/CSU/0_HiRES/data/hires/GSE223917_RAW/`
- **包含**: 7,895 个 `.pairs.gz` 文件

### 已更新的配置文件
1. `stage_files_mapping.csv` - 样本映射文件
2. `20length_chromosome_final_with_path.csv` - 20条染色体数据
3. `21length_chromosome_final_with_path.csv` - 21条染色体数据  
4. `chromosome_counts_final_with_path.csv` - 染色体计数数据
5. `hires_data_processing/contact_tables/20length_chromosome_contact_table.tsv`
6. `hires_data_processing/contact_tables/21length_chromosome_contact_table.tsv`

## Claude 使用注意事项

### 语言设置
**请使用中文回复所有问题和操作说明**

### 常用命令模板
```bash
# 进入项目目录
cd "/Volumes/SumSung500/CSU/0_HiRES"

# 运行 Python 脚本
/Users/wuhaoliu/mamba/envs/schicluster/bin/python your_script.py

# 检查文件
ls -la data/hires/GSE223917_RAW/ | head -10
```

## 项目历史
- 2024-08-30: 完成所有文件路径从旧服务器路径到本地Mac路径的迁移更新
- 路径更新: `/home/duxuyan/Projects/0_HiRES/data/GSE223917_RAW/` → `/Volumes/SumSung500/CSU/0_HiRES/data/hires/GSE223917_RAW/`