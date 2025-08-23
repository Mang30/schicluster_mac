# MPS 加速单细胞 Hi-C 插补使用指南

## 简介

`batch_mps_imputation.py` 提供了使用 Apple Silicon GPU (MPS) 加速的单细胞 Hi-C 数据插补功能。

## 基本用法

### 1. 激活环境
```bash
micromamba activate schicluster
```

### 2. 基本命令
```bash
# 处理单个 stage
python batch_mps_imputation.py --stages E70

# 处理多个 stage  
python batch_mps_imputation.py --stages E70 E75 E80

# 限制处理细胞数量（测试用）
python batch_mps_imputation.py --stages E70 --max-cells 5

# 设置批次大小
python batch_mps_imputation.py --stages E70 --batch-size 5
```

### 3. MPS 相关选项
```bash
# 默认启用 MPS 加速
python batch_mps_imputation.py --stages E70

# 禁用 MPS，使用 CPU
python batch_mps_imputation.py --stages E70 --disable-mps
```

## 参数说明

- `--stages`: 要处理的发育阶段列表 (默认: E70 E75 E80)
- `--max-cells`: 每个 stage 最大处理细胞数 (用于测试)
- `--batch-size`: 批次大小，用于进度管理 (默认: 1)
- `--disable-mps`: 禁用 MPS 加速，使用 CPU
- `--resolution`: Hi-C 数据分辨率 (默认: 20000)

## 输入输出

- **输入**: `/Volumes/SumSung500/CSU/0_HiRES/converted_matrices_by_stage_20k/`
- **输出**: `/Volumes/SumSung500/CSU/0_HiRES/imputed_matrices_by_stage_20k/`
- **日志**: `mps_imputation.log`

## 性能优势

使用 MPS 加速相比 CPU 可以获得 **5-10x** 的速度提升：

- CPU 模式: ~2 秒/染色体
- MPS 模式: ~0.2-0.4 秒/染色体
- 24GB GPU 内存充分利用

## 实际使用示例

### 快速测试
```bash
# 测试单个细胞
python batch_mps_imputation.py --stages E75 --max-cells 1
```

### 批量处理
```bash
# 处理 E70 stage 的所有细胞
python batch_mps_imputation.py --stages E70 --batch-size 10

# 处理所有 stage
python batch_mps_imputation.py --stages E70 E75 E80 E85 E95 EX05 EX15
```

## 监控和日志

- 实时进度显示在控制台
- 详细日志保存在 `mps_imputation.log`
- 自动跳过已处理的文件
- 错误处理和恢复机制

## 故障排除

1. **MPS 不可用**: 检查 PyTorch 版本和 Apple Silicon 支持
2. **内存不足**: 降低批次大小或使用 `--disable-mps`
3. **找不到文件**: 检查输入路径和文件格式

## 输出文件格式

每个染色体生成一个 `.npz` 文件：
```
{cell_name}_{chromosome}_pad1_std1.0_rp0.5_sqrtvc.npz
```

包含插补后的 Hi-C 接触矩阵数据。
