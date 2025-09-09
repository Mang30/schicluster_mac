# JSON to H5AD 转换工具

这是一个将 contact decay profile 数据从 JSON 格式转换为 H5AD 格式的完整工具套件，用于 HiRES 项目的单细胞 Hi-C 数据分析。

## 项目概述

该工具将各个发育阶段（E70, E75, E80, E85, E95, EX05, EX15）的细胞 contact decay profile 数据从 JSON 格式转换为标准的 H5AD 格式，便于后续使用 scanpy 进行单细胞数据分析。

## 功能特性

- ✅ **完整的数据处理流程**: 从 JSON 读取到 H5AD 保存
- ✅ **元数据整合**: 自动匹配和整合 Excel 元数据文件
- ✅ **数据验证**: 全面的数据格式和内容验证
- ✅ **内存优化**: 支持大文件批处理和内存监控
- ✅ **错误处理**: 详细的日志记录和异常处理
- ✅ **进度跟踪**: 实时进度显示和性能监控
- ✅ **配置灵活**: YAML 配置文件支持
- ✅ **批量处理**: 支持单个或多个 stage 批处理

## 文件结构

```
json2h5ad/
├── convert_json_to_h5ad.py    # 核心转换脚本
├── metadata_processor.py     # 元数据处理模块
├── data_validator.py         # 数据验证模块
├── utils.py                  # 工具函数和日志系统
├── config.yaml              # 配置文件
├── process_all_stages.sh     # 批处理脚本
├── README.md                # 使用文档
├── GSE223917_HiRES_emb_metadata.xlsx  # 元数据文件 (需要提供)
├── output/                  # 输出目录
│   ├── stage_E70.h5ad
│   ├── stage_E75.h5ad
│   └── ...
└── logs/                    # 日志目录
    ├── stage_E70_processing_*.log
    └── batch_processing_*.log
```

## 快速开始

### 1. 环境准备

确保已激活正确的 conda 环境：

```bash
# 激活环境
micromamba activate schicluster

# 或直接使用指定的 Python 路径
/Users/wuhaoliu/mamba/envs/schicluster/bin/python
```

### 2. 准备元数据文件

将 `GSE223917_HiRES_emb_metadata.xlsx` 文件放在项目目录中，或在配置文件中指定正确路径。

### 3. 检查配置

编辑 `config.yaml` 文件，确认路径配置：

```yaml
paths:
  input_root: "/Volumes/SumSung500/CSU/0_HiRES/3_create_contact_decay_profile/outputs"
  output_dir: "/Volumes/SumSung500/CSU/0_HiRES/4_contact_decay_profile_2_h5ad/json2h5ad/output"
  metadata_file: "/Volumes/SumSung500/CSU/0_HiRES/4_contact_decay_profile_2_h5ad/json2h5ad/GSE223917_HiRES_emb_metadata.xlsx"
```

## 使用方法

### 方法一：处理单个 Stage

```bash
# 处理单个 stage
/Users/wuhaoliu/mamba/envs/schicluster/bin/python convert_json_to_h5ad.py --stage E75

# 调试模式 (只处理少量文件)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python convert_json_to_h5ad.py --stage E75 --debug

# 指定自定义元数据文件
/Users/wuhaoliu/mamba/envs/schicluster/bin/python convert_json_to_h5ad.py --stage E75 --metadata /Volumes/SumSung500/CSU/0_HiRES/4_contact_decay_profile_2_h5ad/json2h5ad/GSE223917_HiRES_emb_metadata.xlsx
```

### 方法二：批量处理所有 Stages

```bash
# 处理所有 stages
./process_all_stages.sh

# 调试模式
./process_all_stages.sh --debug

# 只处理指定的 stages
./process_all_stages.sh E75 E80 EX05

# 出错时继续处理其他 stages
./process_all_stages.sh --continue
```

### 方法三：使用 make 风格命令 (可选)

创建 `Makefile` 来简化使用：

```makefile
PYTHON = /Users/wuhaoliu/mamba/envs/schicluster/bin/python

# 处理单个 stage
stage-%:
	$(PYTHON) convert_json_to_h5ad.py --stage $*

# 处理所有 stages
all:
	./process_all_stages.sh

# 调试模式
debug:
	./process_all_stages.sh --debug

# 清理输出文件
clean:
	rm -rf output/*.h5ad logs/*.log
```

## 输出说明

### H5AD 文件结构

每个处理后的 H5AD 文件包含：

- **X 矩阵**: 细胞 × 距离 bins 的 contact decay 值矩阵
- **obs (观测值)**: 细胞级别元数据
  - `cell_id`: 细胞标识符
  - `stage`: 发育阶段
  - `celltype`: 细胞类型
  - `cellcycle_phase`: 细胞周期阶段
  - `total_contacts`: 总 contact 数量
  - `power_law_slope`: 幂律衰减斜率
  - 等等...
- **var (变量)**: 特征级别元数据
  - `distance_bins`: 基因组距离 bins
  - `distance_kb`: 以 kb 为单位的距离
- **uns (非结构化)**: stage 级别元数据和处理信息

### 使用 scanpy 读取数据

```python
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# 读取 H5AD 文件
adata = sc.read_h5ad('output/stage_E75.h5ad')

# 查看基本信息
print(f\"细胞数量: {adata.n_obs}\")
print(f\"特征数量: {adata.n_vars}\")

# 查看元数据
print(\"细胞类型分布:\")
print(adata.obs['celltype'].value_counts())

print(\"细胞周期分布:\")
print(adata.obs['cellcycle_phase'].value_counts())

# 可视化 contact decay profile
import numpy as np

# 计算平均 decay profile
mean_decay = np.mean(adata.X, axis=0)
distances_kb = adata.var['distance_kb'].values

plt.figure(figsize=(10, 6))
plt.loglog(distances_kb, mean_decay, 'b-', linewidth=2)
plt.xlabel('Genomic Distance (kb)')
plt.ylabel('Contact Frequency')
plt.title(f'Average Contact Decay Profile - Stage E75')
plt.grid(True, alpha=0.3)
plt.show()

# 按细胞类型分析
for celltype in adata.obs['celltype'].unique():
    if celltype != 'Unknown':
        mask = adata.obs['celltype'] == celltype
        celltype_decay = np.mean(adata.X[mask], axis=0)
        plt.loglog(distances_kb, celltype_decay, label=celltype, linewidth=2)

plt.xlabel('Genomic Distance (kb)')
plt.ylabel('Contact Frequency')
plt.title('Contact Decay Profiles by Cell Type')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

## 配置说明

主要配置选项：

### 路径配置
- `input_root`: 输入数据目录
- `output_dir`: 输出目录
- `metadata_file`: 元数据 Excel 文件路径
- `log_dir`: 日志目录

### 处理参数
- `data_dtype`: 数据类型 (默认 float32)
- `batch_size`: 批处理大小
- `alignment_strategy`: 数据对齐策略

### 验证参数
- `strict_validation`: 是否启用严格验证
- `max_errors`: 允许的最大错误数
- `max_warnings`: 允许的最大警告数

### 输出配置
- `compress_output`: 是否压缩输出文件
- `overwrite_existing`: 是否覆盖已存在的文件
- `create_backup`: 是否创建备份

## 故障排除

### 常见问题

1. **内存不足**
   ```
   解决方案：
   - 减少 batch_size
   - 使用 float32 数据类型
   - 在配置中启用内存优化选项
   ```

2. **元数据匹配失败**
   ```
   检查：
   - JSON 文件名格式是否正确
   - Excel 文件中的 Cellname 字段
   - 细胞 ID 是否一致
   ```

3. **数据验证错误**
   ```
   检查：
   - JSON 文件格式是否正确
   - 必需字段是否存在
   - 数据类型是否正确
   ```

4. **权限问题**
   ```bash
   # 确保脚本有执行权限
   chmod +x process_all_stages.sh
   
   # 确保输出目录可写
   mkdir -p output logs
   chmod 755 output logs
   ```

### 调试技巧

1. **启用调试模式**
   ```bash
   # 只处理少量文件进行测试
   ./process_all_stages.sh --debug
   ```

2. **检查日志文件**
   ```bash
   # 查看最新的日志文件
   ls -lt logs/
   tail -f logs/stage_E75_processing_*.log
   ```

3. **验证环境**
   ```bash
   # 检查 Python 环境
   /Users/wuhaoliu/mamba/envs/schicluster/bin/python -c "import pandas, numpy, scanpy; print('环境正常')"
   ```

## 性能优化

1. **内存优化**
   - 调整 `batch_size` 参数
   - 使用 `float32` 而非 `float64`
   - 启用稀疏矩阵存储 (如果适用)

2. **处理速度**
   - 关闭不必要的验证
   - 减少日志输出级别
   - 使用 SSD 存储

3. **磁盘空间**
   - 启用输出压缩
   - 定期清理日志文件
   - 监控磁盘使用情况

## 贡献和支持

如果遇到问题或需要新功能，请：

1. 检查日志文件中的错误信息
2. 确认配置文件设置正确
3. 验证输入数据格式
4. 查看本文档的故障排除部分

## 版本信息

- **版本**: v1.0
- **创建日期**: 2025-09-09
- **兼容性**: Python 3.8+, scanpy 1.8+
- **测试环境**: macOS, schicluster conda 环境

---

*注意: 请确保在正确的 conda 环境中运行，并且所有必需的 Python 包已安装。*