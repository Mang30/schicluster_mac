# Backup Information

## 原始文件备份

- **原始文件**: `/Volumes/SumSung500/CSU/0_HiRES/schicluster/schicluster/impute/impute_chromosome.py`
- **备份文件**: `/Volumes/SumSung500/CSU/0_HiRES/schicluster/schicluster/impute/impute_chromosome_original.py`
- **备份时间**: 2025年8月22日
- **备份目的**: 在添加 MPS GPU 加速功能之前保护原始代码

## 文件说明

原始文件包含以下主要功能：
1. `calc_sparsity()` - 计算矩阵稀疏度
2. `random_walk_cpu()` - CPU 版本的随机游走算法（性能瓶颈）
3. `impute_chromosome()` - 主要的染色体插补函数

## 待优化部分

- **性能瓶颈**: `random_walk_cpu` 函数（第19-38行）
- **优化方案**: 添加 PyTorch MPS 后端支持，实现 GPU 加速
- **预期提升**: 2-5倍的性能提升

## 恢复方法

如需恢复原始版本，执行：
```bash
cp /Volumes/SumSung500/CSU/0_HiRES/schicluster/schicluster/impute/impute_chromosome_original.py /Volumes/SumSung500/CSU/0_HiRES/schicluster/schicluster/impute/impute_chromosome.py
```
