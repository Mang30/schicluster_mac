# Hi-C接触衰减曲线分析工具

## 项目简介

本项目用于分析Hi-C数据的接触衰减曲线（Contact Decay Profile），支持批量处理多个发育阶段的.cool格式数据，生成接触衰减曲线图表和对比分析结果。

## 功能特性

- 📊 **接触衰减曲线分析**：计算并可视化Hi-C接触频率随基因组距离的衰减
- 🔄 **批量处理**：自动处理多个发育阶段的数据
- ⚡ **并行计算**：支持多进程并行处理，提高分析效率  
- 📈 **多种图表**：生成线性、对数坐标和热图等多种可视化结果
- 🎯 **发育阶段对比**：自动生成各发育阶段间的对比分析
- ⚙️ **灵活配置**：通过YAML配置文件自定义分析参数
- 🐍 **Python生态**：基于cooler、matplotlib等成熟库

## 项目结构

```
3_create_contact_decay_profile/
├── src/                          # 核心代码
│   └── contact_decay_analyzer.py # 接触衰减分析核心类
├── scripts/                      # 脚本工具
│   └── batch_process_stages.py   # 批处理脚本
├── configs/                      # 配置文件
│   └── analysis_config.yaml      # 分析配置
├── utils/                        # 工具模块
│   └── config_loader.py          # 配置加载器
├── outputs/                      # 输出目录
├── logs/                         # 日志目录
├── run_contact_decay_analysis.sh # 主运行脚本
└── README.md                     # 项目说明
```

## 依赖环境

### 必需的Python包
- `numpy >= 1.19.0`
- `pandas >= 1.3.0` 
- `matplotlib >= 3.3.0`
- `seaborn >= 0.11.0`
- `cooler >= 0.8.0`
- `cooltools >= 0.5.0`
- `scipy >= 1.7.0`
- `PyYAML >= 5.4.0`

### Conda环境
推荐使用micromamba的schicluster环境：
```bash
micromamba activate schicluster
```

## 快速开始

### 1. 环境准备
```bash
# 激活conda环境
micromamba activate schicluster

# 安装依赖包（如果尚未安装）
micromamba install -c conda-forge -c bioconda cooler cooltools
pip install seaborn PyYAML
```

### 2. 运行分析
```bash
# 基本运行（每个阶段处理5个文件）
./run_contact_decay_analysis.sh

# 测试模式（每个阶段处理2个文件）
./run_contact_decay_analysis.sh --test

# 处理所有文件
./run_contact_decay_analysis.sh --full

# 自定义参数
./run_contact_decay_analysis.sh -n 10 -w 4 -d 2000
```

### 3. 查看结果
```bash
# 查看批量处理报告
cat outputs/batch_processing_report.json

# 查看生成的图表
ls outputs/*.png

# 查看合并的数据
head outputs/all_stages_decay_profiles.csv
```

## 使用说明

### 单文件分析
```python
from src.contact_decay_analyzer import ContactDecayAnalyzer

# 创建分析器
analyzer = ContactDecayAnalyzer(
    cool_path="path/to/sample.cool",
    output_dir="output_directory"
)

# 运行完整分析
results = analyzer.run_complete_analysis()
```

### 批量处理
```python
from scripts.batch_process_stages import StageWiseDecayProcessor

# 创建批处理器
processor = StageWiseDecayProcessor(
    data_root_dir="../hires_data_processing/outputs",
    output_dir="./outputs"
)

# 运行批量处理
results = processor.run_batch_processing(
    max_files_per_stage=5,
    parallel=True,
    max_workers=4
)
```

### 配置文件使用
```python
from utils.config_loader import get_config

# 加载配置
config = get_config()

# 获取配置值
stages = config.get('stages.stage_list')
max_distance = config.get('analysis.decay_profile.max_distance')
```

## 配置说明

主要配置文件：`configs/analysis_config.yaml`

### 关键配置项
```yaml
# 数据路径
paths:
  data_root: "../hires_data_processing/outputs"
  output_root: "./outputs"

# 发育阶段
stages:
  stage_list: ["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"]

# 分析参数
analysis:
  decay_profile:
    max_distance: 2000  # 最大分析距离(bins)
    use_balanced: true  # 使用平衡化数据

# 并行处理
parallel:
  enable_parallel: true
  max_workers: 4
```

## 输出文件

### 单样本输出
- `{sample}_decay_profile.csv`：衰减曲线数据
- `{sample}_decay_linear.png`：线性坐标图
- `{sample}_decay_loglog.png`：双对数坐标图
- `{sample}_decay_combined.png`：综合对比图
- `{sample}_summary.json`：分析摘要

### 批处理输出
- `batch_processing_report.json`：批处理总报告
- `all_stages_decay_profiles.csv`：所有阶段合并数据
- `stage_comparison_decay_curves.png`：阶段对比曲线图
- `stage_contact_heatmap.png`：接触强度热图
- `stage_{stage}_summary.json`：各阶段摘要

## 命令行参数

### 运行脚本参数
- `-n, --max-files N`：每个阶段最大处理文件数
- `-w, --max-workers N`：并行worker数量
- `-d, --max-distance N`：最大分析距离
- `--full`：处理所有文件
- `--test`：测试模式
- `-h, --help`：显示帮助

### Python脚本参数
```bash
# 单文件分析
python src/contact_decay_analyzer.py -i input.cool -o output_dir

# 批量处理
python scripts/batch_process_stages.py \
    --input data_root \
    --output output_dir \
    --max-files 5 \
    --max-workers 4
```

## 技术原理

### 接触衰减曲线
接触衰减曲线描述了Hi-C接触频率随基因组距离增加而衰减的规律：

1. **距离计算**：计算基因组上任意两点间的线性距离
2. **接触频率**：统计相应距离下的Hi-C接触强度
3. **衰减模式**：拟合幂律衰减模型 P(s) ∝ s^(-α)
4. **可视化**：生成线性和对数坐标下的衰减曲线

### 分析流程
1. 加载.cool格式Hi-C数据
2. 提取接触矩阵（支持染色体特异性分析）
3. 计算不同距离下的平均接触频率
4. 生成衰减曲线和统计图表
5. 进行幂律拟合分析
6. 输出结果文件和可视化图表

## 故障排除

### 常见问题

1. **导入错误**
```bash
# 安装缺失的包
micromamba install -c conda-forge -c bioconda cooler cooltools
pip install seaborn PyYAML
```

2. **内存不足**
```yaml
# 减少并行worker数量
parallel:
  max_workers: 2

# 限制每阶段处理文件数
data_processing:
  max_files_per_stage: 3
```

3. **文件路径错误**
```yaml
# 检查配置文件中的路径
paths:
  data_root: "../hires_data_processing/outputs"  # 相对于项目根目录
```

4. **权限问题**
```bash
# 给脚本添加执行权限
chmod +x run_contact_decay_analysis.sh
```

### 调试模式
```yaml
# 启用调试模式
advanced:
  debugging:
    debug_mode: true
    verbose: true
    save_intermediate: true
```

## 贡献指南

1. Fork项目到个人仓库
2. 创建功能分支：`git checkout -b feature/new-feature`
3. 提交更改：`git commit -am 'Add new feature'`
4. 推送分支：`git push origin feature/new-feature`
5. 创建Pull Request

## 许可证

本项目采用MIT许可证。

## 联系信息

- **作者**：Claude Code Assistant
- **创建时间**：2025-08-29
- **项目类型**：Hi-C数据分析工具

## 更新日志

### v1.0.0 (2025-08-29)
- ✨ 初始版本发布
- 📊 实现接触衰减曲线分析核心功能
- ⚡ 支持批量并行处理
- 📈 提供多种可视化图表
- ⚙️ 完整的配置管理系统
- 📚 详细的文档和使用说明