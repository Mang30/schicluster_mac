# Stage UMAP Plotting Scripts

本目录包含用于绘制各个stage的UMAP图的脚本。

## 文件说明

### 核心脚本

- `plot_stage_umap.py` - 通用的Python UMAP绘图脚本，支持任意stage
- `run_stage_umap.sh` - 通用的Shell启动脚本

### 配置文件
- `color_mapping.json` - 细胞类型颜色映射配置
- `GSE223917_HiRES_emb_metadata.xlsx` - 细胞元数据文件

### 输出目录
- `output/` - 所有stage的UMAP图输出目录
  - `stage_E75/` - E75阶段的UMAP图
  - `stage_E70/` - E70阶段的UMAP图（如果生成）
  - 等等...

## 使用方法

### 方法1: 使用Shell脚本（推荐）

```bash
# 为单个stage绘制UMAP图
./run_stage_umap.sh stage_E75

# 为其他stage绘制UMAP图
./run_stage_umap.sh stage_E70
./run_stage_umap.sh stage_E80

# 指定自定义H5AD文件路径
./run_stage_umap.sh stage_E75 /path/to/custom/stage_E75.h5ad

# 指定自定义输出目录
./run_stage_umap.sh stage_E75 "" /path/to/custom/output
```

### 方法2: 直接使用Python脚本

```bash
# 激活环境
micromamba activate schicluster

# 为单个stage绘制UMAP图
python plot_stage_umap.py --stage stage_E75

# 指定自定义路径
python plot_stage_umap.py --stage stage_E75 --h5ad_file /path/to/file.h5ad --output_dir /path/to/output
```

## 输出文件

每个stage会在`output/<stage_name>/`目录下生成以下文件：

1. `umap_celltype.png` - 按细胞类型着色的UMAP图
2. `umap_stage.png` - 按发育阶段着色的UMAP图
3. `umap_cellcycle.png` - 按细胞周期阶段着色的UMAP图
4. `umap_numeric_features.png` - 按数值特征着色的UMAP图
5. `data_summary.txt` - 数据摘要文件

## 环境要求

- Python 3.12+
- scanpy 1.11+
- pandas 2.3+
- matplotlib
- seaborn
- numpy
- openpyxl

## 输入文件要求

### H5AD文件结构
H5AD文件应包含以下结构：
- `obs` (观测值): 包含细胞元数据
  - `cell_id`: 细胞ID
  - `celltype`: 细胞类型
  - `stage`: 发育阶段
  - `cellcycle_phase`: 细胞周期阶段
  - 其他数值特征...
- `X`: 主数据矩阵（contact decay profile）
- `var`: 变量信息（距离bins）

### 默认输入路径
脚本默认从以下路径读取H5AD文件：
```
/Volumes/SumSung500/CSU/0_HiRES/4_contact_decay_profile_2_h5ad/json2h5ad/output/<stage_name>.h5ad
```

## 颜色映射

细胞类型的颜色通过`color_mapping.json`文件定义，确保不同stage中相同细胞类型使用一致的颜色。

## 错误处理

脚本会自动处理以下问题：
- NaN值：替换为0
- 无限值：替换为有限值
- 数据标准化：根据需要进行
- 维度适配：根据样本数调整PCA和邻居参数

## 示例工作流

1. 确保已生成H5AD文件（通过JSON转换脚本）
2. 运行UMAP绘图脚本：
   ```bash
   ./run_stage_umap.sh stage_E75
   ```
3. 查看输出文件：
   ```bash
   ls -la output/stage_E75/
   ```

## 注意事项

1. 运行前确保conda环境已激活
2. 确保输入的H5AD文件存在且格式正确
3. 大数据集可能需要较长处理时间
4. 确保有足够的内存处理大型数据集
