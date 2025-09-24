# Hi-C Pairs to H5AD 转换器 (内存优化版)

本代码用于将细胞 Hi-C pairs 数据转换为 scanpy 兼容的 h5ad 格式，支持分批处理以优化内存使用。

## 文件说明

- `pairs_to_h5ad_converter.py`: 核心转换器类，支持分批处理
- `run_conversion.py`: 便捷运行脚本，内置内存优化配置
- `README.md`: 使用说明

## 核心特性

### 🧬 雌雄混合细胞支持
- **完整染色体**: 支持全部21条染色体 (chr1-19, chrX, chrY)
- **智能检测**: 自动使用 `mm10_chrom_sizes_with_chrY.txt` 参考文件
- **数据完整性**: 雄性细胞的chrY数据不再丢失
- **一致性处理**: 雌性细胞chrY列为0，雄性细胞包含实际数据

### 🚀 分批处理
- **默认批大小**: 200个细胞/批
- **内存需求**: ~30GB → ~4-8GB
- **处理时间**: 略有增加，但可以在普通机器上运行

### 💾 数据精度优化
- **float32**: 默认选项，节省50%内存
- **float64**: 高精度选项，需要更多内存

### 🎯 接触过滤
- **最小接触数**: 过滤低质量接触，减少数据量
- **默认阈值**: 2次接触

## 使用方法

### 方法1: 便捷脚本 (推荐)

```bash
cd /Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/code

# 默认配置 (推荐：21条染色体，批大小200，float32，过滤≥2次接触)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python run_conversion.py

# 自定义批大小 (适合不同内存配置)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python run_conversion.py --batch_size 100

# 高精度模式
/Users/wuhaoliu/mamba/envs/schicluster/bin/python run_conversion.py --use_float64

# 不分批模式 (需要32-64GB内存)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python run_conversion.py --no_batch

# 使用自动检测染色体模式 (可能丢失chrY数据)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python run_conversion.py --no_chr_file
```

### 方法2: 直接使用转换器

```bash
# 基本用法 (使用染色体参考文件)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python pairs_to_h5ad_converter.py \
    --cell_table ../cell_table.tsv \
    --output ../output/hic_data_100k.h5ad \
    --chr_sizes ../mm10_chrom_sizes_with_chrY.txt \
    --batch_size 200 \
    --min_contacts 2

# 高性能配置 (需要大内存)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python pairs_to_h5ad_converter.py \
    --cell_table ../cell_table.tsv \
    --output ../output/hic_data_100k_highperf.h5ad \
    --chr_sizes ../mm10_chrom_sizes_with_chrY.txt \
    --use_float64 \
    --min_contacts 1 \
    --verbose

# 自动检测模式 (可能不完整)
/Users/wuhaoliu/mamba/envs/schicluster/bin/python pairs_to_h5ad_converter.py \
    --cell_table ../cell_table.tsv \
    --output ../output/hic_data_100k_auto.h5ad \
    --batch_size 200 \
    --verbose
```

```
  # 处理雌性细胞 (20条染色体)
  /Users/wuhaoliu/mamba/envs/schicluster/bin/python
   pairs_to_h5ad_converter.py \
      --cell_table 20length_cell_table.tsv \
      --output ../output/female_cells_20chr.h5ad \
      --batch_size 200

  # 处理雄性细胞 (21条染色体) 
  /Users/wuhaoliu/mamba/envs/schicluster/bin/python
   pairs_to_h5ad_converter.py \
      --cell_table 21length_cell_table.tsv \
      --output ../output/male_cells_21chr.h5ad \
      --batch_size 200
```
## 参数详细说明

### 核心参数
- `--cell_table`: cell_table.tsv 文件路径
- `--output`: 输出 h5ad 文件路径
- `--resolution`: 分辨率（默认: 100000，即100K）
- `--chr_sizes`: 染色体大小文件路径（推荐使用以确保完整性）

### 内存优化参数
- `--batch_size`: 分批处理大小（默认: None，推荐: 100-500）
- `--use_float64`: 使用float64精度（默认: float32）
- `--min_contacts`: 最小接触次数阈值（默认: 1）

### 便捷脚本专用参数
- `--no_chr_file`: 不使用染色体文件，改为自动检测
- `--no_batch`: 禁用分批处理

### 其他参数
- `--verbose`: 详细输出模式

## 内存需求指南

| 配置 | 批大小 | 精度 | 预估内存 | 适用场景 |
|------|--------|------|----------|----------|
| 低内存 | 50 | float32 | 4-8GB | 普通笔记本 |
| 平衡 | 200 | float32 | 8-16GB | 工作站 |
| 高性能 | 500 | float32 | 16-32GB | 高性能服务器 |
| 极速 | 不分批 | float64 | 64-128GB | 大内存服务器 |

## 输出文件命名

输出文件名自动根据配置生成：
```
hic_data_100k_{precision}_{batch_info}_min{min_contacts}_{chr_info}.h5ad
```

示例：
- `hic_data_100k_float32_batch200_min2_ref.h5ad` (使用参考文件，21条染色体)
- `hic_data_100k_float64_nobatch_min1_auto.h5ad` (自动检测，可能20条染色体)
- `hic_data_100k_float32_batch100_min3_ref.h5ad` (完整染色体，低内存配置)

## 输入数据格式

### cell_table.tsv
```
GasaE751001    /path/to/GSM6998595_GasaE751001.pairs.gz
GasaE751002    /path/to/GSM6998596_GasaE751002.pairs.gz
...
```

### pairs.gz 文件格式
标准 Hi-C pairs 格式，包含染色体间相互作用数据：
```
## pairs format v1.0
#sorted: chr1-chr2-pos1-pos2
#chromosome: chr1 195471971
#chromosome: chr2 182113224
...
.    chr1    3052398    chr1    4307873    +    +    .    .
```

## 输出数据

生成的 h5ad 文件包含：
- **X**: 细胞 × 染色体接触特征矩阵 (上三角展平)
- **obs**: 细胞注释信息
- **var**: 特征注释信息 (基因组区间对)
- **uns**: 元数据 (分辨率、染色体信息、处理参数等)

### 元数据 (uns) 包含
- `resolution`: 分辨率
- `chromosomes`: 染色体列表
- `chr_lengths`: 染色体长度
- `bin_mapping`: 分箱映射信息
- `total_bins`: 总分箱数
- `min_contacts`: 接触过滤阈值
- `use_float32`: 数据精度设置

## 性能优化建议

### 内存不足时
1. 减少批大小：`--batch_size 50`
2. 使用float32：默认已启用
3. 提高过滤阈值：`--min_contacts 3`

### 加速处理时
1. 增加批大小：`--batch_size 500`
2. 使用SSD存储
3. 关闭其他程序释放内存

## 依赖包

确保 schicluster 环境中安装了以下包：
- pandas
- numpy
- scanpy
- anndata
- scipy
- tqdm

## 故障排除

### 常见错误
1. **内存不足**: 减少batch_size或使用float32
2. **文件路径错误**: 检查cell_table.tsv中的路径
3. **磁盘空间不足**: 确保输出目录有足够空间

### 监控进度
- 使用 `--verbose` 参数查看详细日志
- 每100个细胞会输出进度信息
- 分批模式下每批完成后显示进度

## 注意事项

1. **内存规划**: 根据可用内存选择合适的batch_size
2. **时间预估**: 7,469个细胞预计需要2-6小时
3. **存储空间**: 输出文件可能达到10-50GB
4. **稳定运行**: 建议在稳定的服务器环境中运行
5. **后续分析**: h5ad文件可直接用于scanpy分析