# 🚀 一键启动脚本使用说明

## 📋 概述

这个一键启动脚本 (`run.sh`) 提供了完整的 Hi-C 数据处理流水线，支持：
- 自动发现和处理所有可用的 stages
- 智能硬件检测 (MPS/CUDA/CPU)
- 完整流水线：矩阵合并 → h5ad构建 → UMAP可视化
- 灵活的参数配置和跳步控制

## 🎯 快速开始

### 1. 基本功能测试
```bash
# 快速验证功能（推荐首次使用）
./run.sh --mode test --use-distance

# 仅测试 h5ad 构建（跳过其他步骤）
./run.sh --mode test --test-cells 2 --use-distance --skip-merge --skip-umap
```

### 2. 处理指定 stages
```bash
# 处理特定的 stages
./run.sh --mode stage --stages E70,E80 --use-distance

# 限制每个 stage 的细胞数
./run.sh --mode stage --stages E70 --max-cells 10 --use-distance
```

### 3. 全量处理（生产环境）
```bash
# 处理所有可用的 stages（推荐用于正式分析）
./run.sh --mode all --use-distance

# 强制使用 CPU（如果GPU有问题）
./run.sh --mode all --use-distance --gpu-type cpu
```

## 📊 运行模式详解

### 🧪 测试模式 (`--mode test`)
- 用途：快速功能验证
- 处理范围：前3个发现的stages
- 默认细胞数：每stage 3个细胞
- 推荐场景：首次使用、调试、功能验证

### 🎯 指定模式 (`--mode stage`)
- 用途：处理特定的stages
- 处理范围：`--stages` 参数指定
- 细胞数：可通过 `--max-cells` 限制
- 推荐场景：分析特定发育阶段、分步处理

### 🚀 全量模式 (`--mode all`)
- 用途：生产环境批量处理
- 处理范围：所有发现的stages
- 细胞数：默认处理全部细胞
- 推荐场景：正式分析、论文数据生成

## ⚡ 性能优化建议

### 1. 使用距离特征（强烈推荐）
```bash
--use-distance    # 4K特征 vs 114M特征，速度提升5-10倍
```

### 2. GPU加速选择
```bash
--gpu-type auto   # 自动检测（默认）
--gpu-type mps    # Apple Silicon (推荐Mac用户)
--gpu-type cuda   # NVIDIA GPU
--gpu-type cpu    # CPU模式（兼容性最好）
```

### 3. 分步处理
```bash
# 先合并矩阵
./run.sh --mode all --skip-h5ad --skip-umap

# 再构建h5ad
./run.sh --mode all --skip-merge --skip-umap --use-distance

# 最后生成UMAP
./run.sh --mode all --skip-merge --skip-h5ad
```

## 📁 输出文件组织

```
output/
├── merged_matrices/          # 合并后的细胞矩阵
│   ├── E70/
│   │   ├── GasfE703172_merged.npz
│   │   └── ...
│   └── ...
├── hic_h5ad_files/          # 处理后的h5ad文件
│   ├── E70_processed.h5ad
│   ├── E80_processed.h5ad
│   └── ...
├── umap_plots/              # UMAP可视化图表
│   ├── E70_umap.png
│   ├── E80_umap.png
│   └── ...
├── logs/                    # 运行日志
│   └── run.log
└── reports/                 # 处理报告
    └── summary.json
```

## 🔧 高级用法

### 自定义处理流程
```bash
# 仅重新生成UMAP（h5ad文件已存在）
./run.sh --mode all --skip-merge --skip-h5ad

# 从矩阵合并开始，但不生成UMAP
./run.sh --mode all --skip-umap --use-distance

# 仅处理特定stage的UMAP
./run.sh --mode stage --stages E70 --skip-merge --skip-h5ad
```

### 批量处理不同参数组合
```bash
# 脚本1: 快速预览（距离特征）
./run.sh --mode test --use-distance

# 脚本2: 完整特征分析
./run.sh --mode stage --stages E70 --max-cells 5

# 脚本3: 生产环境处理
./run.sh --mode all --use-distance
```

## 📊 性能基准

### 处理时间估算
- **测试模式**: 3 stages × 3 cells ≈ 2-5分钟
- **单stage处理**: 10 cells ≈ 1-2分钟 (距离特征)
- **全量处理**: 取决于数据规模，通常 10-60分钟

### 硬件要求
- **内存**: 建议 8GB+ (距离特征) / 32GB+ (完整特征)
- **存储**: 每个stage约 100MB-1GB 输出文件
- **GPU**: MPS(Apple Silicon) > CUDA > CPU

## 🐛 常见问题

### 1. 路径问题
- 确保在 `plotumap` 目录下运行脚本
- 检查 `config/` 目录是否包含配置文件
- 验证 `core/` 目录包含所有Python模块

### 2. 内存不足
```bash
# 使用距离特征减少内存使用
./run.sh --mode test --use-distance

# 分批处理
./run.sh --mode stage --stages E70 --max-cells 5
```

### 3. GPU问题
```bash
# 强制使用CPU
./run.sh --mode test --gpu-type cpu

# 检查GPU可用性
python3 -c "import torch; print(f'MPS: {torch.backends.mps.is_available()}')"
```

### 4. Python依赖问题
```bash
# 确保在正确的conda环境中
conda activate schicluster

# 检查关键包
python3 -c "import torch, scanpy, anndata; print('Dependencies OK')"
```

## 📝 日志和监控

### 查看运行日志
```bash
# 实时查看日志
tail -f output/logs/run.log

# 查看处理总结
cat output/reports/summary.json
```

### 进度监控
脚本提供实时进度条显示：
```
[进度] [E70] 构建h5ad: [████████████░░░░░░░░] 67% (2/3)
```

---

**提示**: 建议首次使用时运行 `./run.sh --mode test --use-distance` 进行功能验证！ 🎯
