# Contact Decay Profile JSON to H5AD Converter

## 项目概述

4_contact_decay_profile_2_h5ad/json2h5ad项目将处理后的各个stage下的细胞contact decay profile数据从JSON格式转换为H5AD格式，便于后续的单细胞数据分析。

## 输入数据规范

### 数据路径

```text
3_create_contact_decay_profile/outputs/
├── stage_E70/
├── stage_E75/
├── stage_E80/
├── stage_E85/
├── stage_E95/
├── stage_EX05/
└── stage_EX15/
```

### 元数据文件

- **文件路径**: `4_contact_decay_profile_2_h5ad/json2h5ad/GSE223917_HiRES_emb_metadata.xlsx`
- **用途**: 提供细胞类型、发育阶段、细胞周期等重要元数据信息
- **关键字段**:
  - `Cellname`: 细胞名称/ID (如 "GasaE751001")
  - `Stage`: 发育阶段 (E70, E75, E80, E85, E95, EX05, EX15)
  - `Celltype`: 细胞类型 (如 "ExE ectoderm", "neural ectoderm", "ExE mesoderm", "early mesenchyme"等)
  - `Cellcycle phase`: 细胞周期阶段 (G0, G1, Early-S, Mid-S, Late-S, G2, M)
  - `CDPS cluster`: CDPS聚类标签
  - `Sub_k_cluster`: 子聚类标签
  - `Raw contacts`, `Clean3 contacts`: 接触数量统计
  - `G1S.Score`, `G2M.Score`: 细胞周期评分
  - `repli score`: 复制评分

### 输入文件格式

- **文件类型**: JSON格式
- **文件命名**: 按细胞ID命名的contact decay profile文件（例如：`GasdE701023_decay_profile.json`）
- **主要数据文件结构**:

```json
{
  "distances": [2, 3, 4, 5, 6, ...],              // 基因组距离（bins）
  "decay_values": [0.0314, 0.0240, 0.0173, ...],  // Contact decay值
  "distance_kb": [200.0, 300.0, 400.0, ...],      // 以kb为单位的距离
  "total_contacts": 10,                            // 总contact数
  "resolution": 100000,                            // 分辨率(bp)
  "power_law_slope": -1.0452775367553244           // 幂律斜率
}
```

- **摘要文件结构**（`analysis_summary.json`）:

```json
{
  "sample_name": "GasfE703185",
  "cool_file": "/path/to/file.cool",
  "resolution": 100000,
  "total_bins": 26350,
  "power_law_slope": -0.9296874402856548,
  "data_points": 24713
}
```

## 输出数据规范

### 输出路径

```text
4_contact_decay_profile_2_h5ad/json2h5ad/output/
├── stage_E70.h5ad
├── stage_E75.h5ad
├── stage_E80.h5ad
├── stage_E85.h5ad
├── stage_E95.h5ad
├── stage_EX05.h5ad
└── stage_EX15.h5ad
```

### H5AD文件结构要求

- **观测值 (obs)**: 细胞级别的元数据信息
  - `cell_id`: 细胞标识符（如 "GasdE701023"）
  - `stage`: 发育阶段（如 "E70"）
  - `celltype`: 细胞类型（如 "ExE ectoderm", "neural ectoderm"等）
  - `cellcycle_phase`: 细胞周期阶段（如 "G1", "S", "G2", "M"）
  - `cdps_cluster`: CDPS聚类标签
  - `sub_k_cluster`: 子聚类标签
  - `total_contacts`: 总contact数量
  - `power_law_slope`: 幂律衰减斜率
  - `resolution`: 数据分辨率
  - `g1s_score`: G1/S期评分（如果可用）
  - `g2m_score`: G2/M期评分（如果可用）
  - `repli_score`: 复制评分（如果可用）
  
- **变量 (var)**: 距离bins特征
  - `distance_bins`: 基因组距离bins
  - `distance_kb`: 以kb为单位的距离
  
- **主数据矩阵 (X)**: 细胞 × 距离bins的contact decay值矩阵

- **非结构化数据 (uns)**: stage级别的元数据和处理信息
  - `processing_date`: 处理日期
  - `stage_summary`: 每个stage的细胞数量统计
  - `celltype_counts`: 细胞类型分布统计

## 数据处理特殊要求

### 文件处理规则

1. **过滤文件**: 忽略 `analysis_summary.json` 文件，只处理 `*_decay_profile.json` 文件
2. **细胞ID提取**: 从文件名中提取细胞ID（去除 `_decay_profile.json` 后缀）
3. **元数据匹配**: 根据细胞ID从Excel元数据文件中匹配细胞类型和stage信息
4. **数据对齐**: 确保所有细胞的距离bins一致，处理长度不同的情况
5. **缺失值处理**: 对于不同长度的decay_values，使用NaN填充或截断到最小长度

### 元数据处理规则

1. **ID匹配策略**: 
   - JSON文件名格式: `GasdE701023_decay_profile.json` -> 细胞ID: `GasdE701023`
   - Excel中的Cellname格式: `GasaE751001`, `GasdE701023` 等
   - 处理可能的ID格式不一致问题
2. **缺失元数据处理**: 
   - 如果细胞在Excel中找不到，使用默认值或标记为"Unknown"
   - 保留原始文件名作为备用标识
3. **数据类型转换**: 
   - 确保stage信息格式一致（如"E70"）
   - 细胞类型名称标准化
   - 数值型元数据的适当处理

### 数据验证规则

1. **必需字段**: 确保每个JSON文件包含 `distances`, `decay_values`, `distance_kb` 字段
2. **数据一致性**: 验证 `distances`, `decay_values`, `distance_kb` 长度一致
3. **数值范围**: 检查 `decay_values` 为非负数，`distances` 为递增序列
4. **元数据完整性**: 验证关键元数据字段的存在和有效性

### 内存优化

- 分批处理大量JSON文件，避免一次性加载所有数据
- 使用稀疏矩阵存储如果decay_values中有大量零值
- 合理设置数据类型（float32 vs float64）
- 预先加载并缓存元数据Excel文件

## 技术要求

### 依赖库

```python
import pandas as pd
import numpy as np
import scanpy as sc
import json
import os
from pathlib import Path
import logging
import openpyxl  # 用于读取Excel文件
from datetime import datetime
import warnings
```

### 环境配置

- **Conda环境**: `micromamba activate schicluster`
- **Python路径**: `/Users/wuhaoliu/mamba/envs/schicluster/bin/python`

## 代码结构要求

### 核心模块

1. **元数据读取模块**: 读取和解析Excel元数据文件
2. **数据读取模块**: 读取JSON文件并解析数据
3. **数据匹配模块**: 将JSON数据与元数据进行匹配整合
4. **数据转换模块**: 将整合后的数据转换为AnnData格式
5. **质量控制模块**: 数据验证和清洗
6. **输出模块**: 保存H5AD文件

### 脚本要求

1. **核心转换脚本**: `convert_json_to_h5ad.py`
2. **批处理脚本**: 7个stage对应的启动脚本
3. **配置文件**: 存储路径和参数配置
4. **元数据处理脚本**: `process_metadata.py`（可选独立模块）

## 错误处理要求

### 日志记录

- 使用Python logging模块
- 记录处理进度和错误信息
- 日志文件按stage分别保存

### 异常处理

- 文件不存在或损坏的处理
- 数据格式不符合预期的处理
- 内存不足的处理

## 质量控制

### 数据验证

- 检查输入文件完整性
- 验证数据格式正确性
- 确保输出文件结构符合规范

### 性能要求

- 支持大文件处理
- 内存优化
- 进度显示

## 使用示例

### 单个stage处理

```bash
/Users/wuhaoliu/mamba/envs/schicluster/bin/python convert_json_to_h5ad.py --stage E75 --metadata GSE223917_HiRES_emb_metadata.xlsx
```

### 批量处理

```bash
bash process_all_stages.sh
```

### 数据读取示例

```python
import scanpy as sc
import pandas as pd

# 读取生成的H5AD文件
adata = sc.read_h5ad('output/stage_E75.h5ad')

# 查看元数据
print("细胞数量:", adata.n_obs)
print("特征数量:", adata.n_vars)
print("细胞类型分布:")
print(adata.obs['celltype'].value_counts())
print("发育阶段:")
print(adata.obs['stage'].unique())
print("细胞周期分布:")
print(adata.obs['cellcycle_phase'].value_counts())
```

## 文件命名规范

- 脚本文件: `snake_case.py`
- 配置文件: `config.yaml`
- 日志文件: `stage_name_processing_YYYYMMDD_HHMMSS.log`

## 版本信息

- **创建日期**: 2025-09-09
- **版本**: v1.0
- **最后更新**: 2025-09-09
- **维护者**: Claude AI Assistant

## 注意事项

1. **元数据匹配**: 确保JSON文件名中的细胞ID能正确匹配Excel文件中的Cellname
2. **数据完整性**: 处理前确认输入数据和元数据文件的完整性
3. **内存管理**: 注意内存使用，大文件分批处理
4. **数据一致性**: 保持输出文件的一致性格式
5. **版本控制**: 定期备份重要数据和元数据文件
6. **缺失值策略**: 明确定义如何处理缺失的元数据信息
7. **数据验证**: 转换完成后验证H5AD文件的元数据完整性和正确性

## 元数据字段映射

| Excel字段 | H5AD obs字段 | 说明 |
|-----------|--------------|------|
| Cellname | cell_id | 细胞标识符 |
| Stage | stage | 发育阶段 |
| Celltype | celltype | 细胞类型 |
| Cellcycle phase | cellcycle_phase | 细胞周期阶段 |
| CDPS cluster | cdps_cluster | CDPS聚类标签 |
| Sub_k_cluster | sub_k_cluster | 子聚类标签 |
| G1S.Score | g1s_score | G1/S期评分 |
| G2M.Score | g2m_score | G2/M期评分 |
| repli score | repli_score | 复制评分 |
| Raw contacts | raw_contacts | 原始接触数 |
| Clean3 contacts | clean_contacts | 清洗后接触数 |
