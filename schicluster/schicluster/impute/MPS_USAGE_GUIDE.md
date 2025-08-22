# MPS GPU 加速插补使用指南

## 概述
修改后的 `impute_chromosome.py` 现在支持 Apple Silicon Mac 上的 MPS (Metal Performance Shaders) GPU 加速。

## 新增功能

### 1. 自动检测和回退
- 自动检测是否安装了 PyTorch 和 MPS 支持
- 如果 MPS 不可用，自动回退到 CPU 实现
- 详细的日志记录显示使用的加速方法

### 2. MPS 加速的组件
- **Random Walk with Restart**: 使用 GPU 稀疏矩阵运算
- **Gaussian Convolution**: 使用 GPU 卷积操作

### 3. 新增参数
```python
def impute_chromosome(..., use_mps=True):
```
- `use_mps`: 是否启用 MPS 加速（默认 True）

## 安装要求

### 安装 PyTorch with MPS 支持
```bash
pip install torch torchvision torchaudio
```

### 验证 MPS 可用性
```python
import torch
print("MPS available:", torch.backends.mps.is_available())
```

## 使用示例

### 基本使用（自动GPU加速）
```python
from schicluster.impute.impute_chromosome import impute_chromosome

# 默认启用 MPS 加速
impute_chromosome(
    chrom='chr1',
    resolution=20000,
    output_path='output.npz',
    contact_path='contacts.txt',
    chrom_size_path='chrom_sizes.txt'
)
```

### 强制使用 CPU
```python
# 禁用 MPS，强制使用 CPU
impute_chromosome(
    chrom='chr1',
    resolution=20000,
    output_path='output.npz',
    contact_path='contacts.txt',
    chrom_size_path='chrom_sizes.txt',
    use_mps=False
)
```

## 性能优势

### 预期加速比
- **Random Walk**: 2-5倍加速（取决于矩阵大小）
- **Gaussian Filter**: 1.5-3倍加速
- **整体流程**: 1.5-3倍加速

### 最佳性能条件
- 大矩阵（> 10000 x 10000）
- 高迭代次数的 Random Walk
- 大的 Gaussian 卷积核

## 日志输出示例

```
INFO: MPS acceleration is available and will be used for GPU acceleration.
INFO: Starting chromosome chr1 imputation with MPS GPU acceleration
INFO: Using MPS-accelerated Gaussian convolution
INFO: Using MPS acceleration for Random Walk with Restart
DEBUG: MPS Iter 1 takes 0.123 seconds. Loss: 0.045; Sparsity: 0.002
```

## 故障排除

### 常见问题

1. **PyTorch 未安装**
   ```
   INFO: PyTorch not found. Using CPU-only implementation.
   ```
   解决方案：安装 PyTorch

2. **MPS 不可用**
   ```
   INFO: MPS acceleration is not available. Falling back to CPU.
   ```
   解决方案：确保使用 Apple Silicon Mac 和支持 MPS 的 PyTorch 版本

3. **内存不足**
   ```
   WARNING: MPS acceleration failed: ... Falling back to CPU.
   ```
   解决方案：使用 `use_mps=False` 或减小矩阵大小

### 调试技巧
- 启用详细日志：`logging.basicConfig(level=logging.DEBUG)`
- 监控 GPU 内存使用：Activity Monitor > GPU History
- 比较 CPU vs MPS 性能：运行相同任务测试时间

## 兼容性

### 支持的系统
- ✅ Apple Silicon Mac (M1/M2/M3)
- ✅ Intel Mac (CPU only)
- ✅ Linux (CPU only)
- ✅ Windows (CPU only)

### 向后兼容性
- 完全向后兼容原有 API
- 原有调用代码无需修改
- 新参数有合理默认值
