# 使用说明：为h5ad文件添加元数据并绘制UMAP图

## 任务概述

### 任务1: 为h5ad文件添加元数据 (1_add_metadata_to_h5ad.py)
- **输入**: 4_contact_decay_profile_2_h5ad/outputs目录中的h5ad文件
- **元数据来源**: 5_plot_umap/GSE223917_HiRES_emb_metadata.xlsx
- **输出**: 4_contact_decay_profile_2_h5ad/outputs_with_metadata目录中更新后的h5ad文件

### 任务2: 绘制UMAP图 (2_plot_umap_from_h5ad.py)
- **输入**: 带元数据的h5ad文件
- **颜色映射**: 5_plot_umap/color_mapping.json
- **输出**: 5_plot_umap/umap_plots目录中的PNG图像文件

## 执行顺序

1. 首先运行任务1脚本添加元数据：
   ```bash
   cd /home/duxuyan/Projects/schicluster_mac/5_plot_umap
   ./1_run_add_metadata.sh
   ```

2. 然后运行任务2脚本绘制UMAP图：
   ```bash
   cd /home/duxuyan/Projects/schicluster_mac/5_plot_umap
   ./2_run_plot_umap.sh
   ```

## 文件说明

### 任务1相关文件
- `1_add_metadata_to_h5ad.py`: Python脚本，负责将元数据添加到h5ad文件
- `1_run_add_metadata.sh`: Bash脚本，用于执行Python脚本

### 任务2相关文件
- `2_plot_umap_from_h5ad.py`: Python脚本，负责从h5ad文件绘制UMAP图
- `2_run_plot_umap.sh`: Bash脚本，用于执行Python脚本

## 输出结果

### 任务1输出
- `4_contact_decay_profile_2_h5ad/outputs_with_metadata/stage_E75_decay_profiles.h5ad`
- `4_contact_decay_profile_2_h5ad/outputs_with_metadata/stage_E80_decay_profiles.h5ad`

### 任务2输出
- `5_plot_umap/umap_plots/stage_E75_decay_profiles_umap_celltype.png`
- `5_plot_umap/umap_plots/stage_E75_decay_profiles_umap_stage.png`
- `5_plot_umap/umap_plots/stage_E80_decay_profiles_umap_celltype.png`
- `5_plot_umap/umap_plots/stage_E80_decay_profiles_umap_stage.png`

## 注意事项

1. 确保已安装micromamba并配置好`3_schicluster_python38`环境
2. 脚本会自动检查并安装必要的Python包
3. 颜色映射确保相同细胞类型在不同stage中使用相同颜色
4. 如果h5ad文件中没有预计算的UMAP坐标，脚本会自动计算