# Hi-C NPZ 到 AnnData 转换工具

这个工具包提供了将单细胞 Hi-C 数据从 NPZ 格式转换为 AnnData 格式的完整流水线。

## 功能特点

- **模块化设计**: 各个组件可独立使用和测试
- **高效处理**: 支持稀疏矩阵的高效处理
- **批量转换**: 支持按发育阶段批量处理数据
- **灵活配置**: 可自定义染色体、文件模式等参数
- **完整日志**: 详细的处理日志和错误报告

## 文件结构

```
npz2h5ad/
├── hic_converter.py          # 核心工具库
├── convert_npz_to_h5ad.py    # 主转换程序
├── test_conversion.py        # 测试脚本
├── README.md                 # 说明文档
└── output/                   # 输出目录（自动创建）
```

## 核心组件

### 1. HiCNpzLoader
负责从 NPZ 文件加载稀疏矩阵，支持 COO 和 CSR 格式。

### 2. UpperTriangleExtractor
提取矩阵的上三角部分并展平为一维向量，可选择是否包含对角线元素。

### 3. ChromosomeManager
管理染色体顺序和文件验证，确保数据的一致性。

### 4. FeatureGenerator
为单个细胞生成特征向量，将所有染色体的上三角矩阵拼接。

### 5. AnnDataBuilder
构建 AnnData 对象，包含完整的元数据信息。

## 使用方法

### 1. 测试工具功能

首先运行测试脚本确保工具正常工作：

```bash
cd npz2h5ad
python test_conversion.py
```

### 2. 转换所有阶段数据

```bash
python convert_npz_to_h5ad.py \
    --input_dir ../output/imputed_matrices_by_stage \
    --output_dir ./output
```

### 3. 转换特定阶段

```bash
python convert_npz_to_h5ad.py \
    --input_dir ../output/imputed_matrices_by_stage \
    --output_dir ./output \
    --stages E70,E75,E80
```

### 4. 自定义参数

```bash
python convert_npz_to_h5ad.py \
    --input_dir ../output/imputed_matrices_by_stage \
    --output_dir ./output \
    --no_diagonal \
    --chromosomes chr1,chr2,chr3 \
    --max_cells_per_stage 50 \
    --verbose
```

## 参数说明

- `--input_dir`: 输入目录，包含按阶段组织的细胞数据
- `--output_dir`: 输出目录，保存 H5AD 文件
- `--stages`: 要处理的发育阶段，逗号分隔
- `--include_diagonal`/`--no_diagonal`: 是否包含对角线元素
- `--file_pattern`: 染色体文件的命名模式
- `--max_cells_per_stage`: 每个阶段处理的最大细胞数（用于测试）
- `--chromosomes`: 要包含的染色体列表
- `--verbose`: 启用详细日志

## 输出格式

每个发育阶段生成一个 H5AD 文件，包含：

### 观测数据 (adata.obs)
- `cell_id`: 细胞标识符
- `stage`: 发育阶段

### 特征数据 (adata.var)
- `chromosome`: 染色体名称
- `bin1`: 第一个 bin 的索引
- `bin2`: 第二个 bin 的索引
- `feature_name`: 特征名称 (chr_binI_binJ 格式)

### 表达矩阵 (adata.X)
- 形状: (n_cells, n_features)
- 内容: 每个细胞所有染色体上三角矩阵的拼接向量

### 元数据 (adata.uns)
- `processing_info`: 处理信息
  - `conversion_method`: 转换方法
  - `include_diagonal`: 是否包含对角线
  - `n_chromosomes`: 染色体数量
  - `chromosomes`: 染色体列表

## 内存优化

- 对于大矩阵（>100k 非零元素），自动使用稀疏矩阵优化算法
- 按细胞逐个处理，避免内存溢出
- 支持限制细胞数量进行测试

## 错误处理

- 详细的日志记录和错误报告
- 自动跳过缺失的染色体文件
- 继续处理即使某些细胞失败

## 注意事项

1. 确保输入数据按照预期的目录结构组织
2. 所有 NPZ 文件应包含方形矩阵（用于上三角提取）
3. 不同细胞的矩阵维度应该一致
4. 建议先用小数据集测试

## 依赖要求

```
numpy
scipy
pandas
anndata
scanpy (可选，用于后续分析)
```

## 示例输出

```
INFO - Processing stage: E70
INFO - Found 384 cells in stage E70
INFO - Processing cell 1/384: GasdE701001
INFO - Processed chr1: shape (2000, 2000), features 2001000
INFO - Processed chr2: shape (1800, 1800), features 1621800
...
INFO - Created AnnData object: AnnData object with n_obs × n_vars = 384 × 45678900
INFO - Stage E70 summary:
INFO -   Cells: 384
INFO -   Features: 45678900
INFO -   Chromosomes: 20
INFO -   Output file: ./output/E70_hic_data.h5ad
```
