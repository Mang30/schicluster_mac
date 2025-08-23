# scHiCluster HiC数据插补指南

## 概述

本指南介绍如何使用schicluster对HiC数据进行插补，数据按发育阶段（Stage）组织。

## 数据分析结果

你的GSE223917_RAW数据包含：
- **总文件数**: 7,895个pairs.gz文件
- **成功识别的文件**: 7,469个（94.6%）
- **发育阶段**（共7个Stage）:
  - **E70**: 557个细胞
  - **E75**: 1,870个细胞  
  - **E80**: 559个细胞
  - **E85**: 1,125个细胞
  - **E95**: 1,108个细胞
  - **EX05**: 1,128个细胞
  - **EX15**: 1,122个细胞

## 文件结构

```
lwh_code/
├── convert_pairs_to_matrix.py          # 基础转换脚本
├── stage_aware_conversion.py           # 按Stage组织的转换脚本
├── stage_aware_imputation.sh          # 按Stage组织的插补脚本
├── batch_convert_pairs.py              # 批量转换脚本
└── README_schicluster_imputation.md    # 本文档

converted_matrices_by_stage/            # 转换后的稀疏矩阵（按Stage组织）
├── E70/          # 557个细胞
├── E75/          # 1,870个细胞  
├── E80/          # 559个细胞
├── E85/          # 1,125个细胞
├── E95/          # 1,108个细胞
├── EX05/         # 1,128个细胞
└── EX15/         # 1,122个细胞

imputed_matrices_by_stage/              # 插补后的数据（按Stage组织）
├── E70/
├── E75/
├── E80/
├── E85/
├── E95/
├── EX05/
└── EX15/
```

## 使用步骤

### 1. 激活schicluster环境

```bash
micromamba activate schicluster
```

### 2. 数据转换（按Stage组织）

转换少量数据测试：
```bash
python lwh_code/stage_aware_conversion.py \
    --input_dir "data/GSE223917_RAW" \
    --output_dir "converted_matrices_by_stage" \
    --chrom_sizes "mm10_chrom_sizes.txt" \
    --max_per_stage 5
```

转换所有数据：
```bash
python lwh_code/stage_aware_conversion.py \
    --input_dir "data/GSE223917_RAW" \
    --output_dir "converted_matrices_by_stage" \
    --chrom_sizes "mm10_chrom_sizes.txt"
```

### 3. 运行插补

确保schicluster环境已激活，然后运行：

```bash
chmod +x lwh_code/stage_aware_imputation.sh
bash lwh_code/stage_aware_imputation.sh
```

## 参数说明

### 转换参数
- `--resolution`: 分辨率，默认100kb
- `--max_per_stage`: 每个Stage最大处理文件数（用于测试）

### 插补参数
- `PAD=1`: 高斯卷积的截断参数
- `STD=1`: 高斯卷积的标准差
- `RP=0.5`: 随机游走重启概率

## 输出文件

### 转换后文件格式
每个细胞每个染色体生成一个稀疏矩阵文件：
```
{cell_id}_{chromosome}.txt
```
内容格式：`bin1	bin2	count`

### 插补后文件格式
schicluster会生成：
- `.npz`文件：插补后的稀疏矩阵
- 其他中间文件

## 数据组织优势

按Stage组织数据的优势：
1. **便于分析**: 可以比较不同发育阶段的3D基因组结构
2. **并行处理**: 可以针对不同Stage并行处理  
3. **内存友好**: 避免一次性加载所有数据
4. **结果清晰**: 插补结果按发育阶段组织，便于后续分析
5. **高覆盖率**: 94.6%的数据文件被正确识别和分类

## 快速数据概览

运行以下命令查看完整数据统计：
```bash
python lwh_code/data_statistics.py
```

## 下一步分析

插补完成后，可以进行：
1. **聚类分析**: 使用插补后的数据进行细胞聚类
2. **差异分析**: 比较不同Stage之间的3D基因组差异
3. **可视化**: 使用schicluster的可视化功能
4. **特征提取**: 计算compartment、domain、loop等特征

## 注意事项

1. **内存需求**: 插补过程需要大量内存，建议在高配置服务器上运行
2. **时间消耗**: 全量数据插补可能需要数小时到数天
3. **磁盘空间**: 确保有足够的磁盘空间存储插补结果
4. **环境依赖**: 确保schicluster环境正确安装和激活

## 故障排除

如果遇到问题：
1. 检查schicluster环境是否正确激活
2. 检查输入文件路径是否正确
3. 检查磁盘空间是否充足
4. 查看错误日志定位问题

## 文件说明

- `mm10_chrom_sizes.txt`: 小鼠基因组染色体大小文件
- 脚本自动从pairs文件头部提取染色体信息生成此文件