# Hi-C 单细胞数据 h5ad 构建 - 性能对比与使用指南

## 📊 性能测试结果

### 测试环境
- **系统**: macOS (Apple Silicon)
- **测试数据**: E70 stage，3个细胞，20条染色体
- **硬件**: MPS GPU 支持

### 方法对比

| 方法 | 特征维度 | 处理3细胞时间 | 预估400细胞时间 | 内存使用 | 推荐度 |
|------|----------|---------------|----------------|----------|--------|
| **CPU上三角** | 1.15亿维 | >5分钟 | 8-20小时 | >10GB | ❌ 不推荐 |
| **MPS距离特征** | 4000维 | 1分32秒 | 3-4小时 | <0.1GB | ✅ **强烈推荐** |
| **CUDA加速** | 1.15亿维 | 预计2-3分钟 | 4-6小时 | 2-5GB | ⚠️ 需NVIDIA GPU |

## 🚀 推荐使用方案

### 方案1: MPS距离特征 (推荐)
```bash
python build_stage_h5ad_mps.py \
  --stage-dir /path/to/E70 \
  --obs-xlsx /path/to/metadata.xlsx \
  --output /path/to/E70_dist.h5ad \
  --batch-size 8 \
  --distance-features
```

**优点**:
- ⚡ 速度快 (1.5分钟/3细胞)
- 💾 内存友好 (<0.1GB)
- 🎯 保留Hi-C关键信息
- 📱 适用于Apple Silicon

**特征内容**: 每条染色体的距离分箱统计 (mean, std, sum, count)

### 方案2: CUDA完整上三角 (高性能机器)
```bash
# 需要先安装: pip install cupy-cuda11x  # 或对应CUDA版本
python build_stage_h5ad_cuda.py \
  --stage-dir /path/to/E70 \
  --obs-xlsx /path/to/metadata.xlsx \
  --output /path/to/E70_full.h5ad \
  --batch-size 16
```

**适用场景**:
- 🖥️ 有NVIDIA GPU
- 🔬 需要完整上三角信息
- 💪 大内存机器 (16GB+)

### 方案3: CPU简化版 (备选)
```bash
python build_stage_h5ad_simple.py \
  --stage-dir /path/to/E70 \
  --obs-xlsx /path/to/metadata.xlsx \
  --output /path/to/E70_simple.h5ad \
  --max-cells 50  # 限制细胞数
```

**适用场景**:
- 🕐 不急于处理
- 📊 只需要少量细胞做测试
- 💻 普通CPU机器

## 📈 批量处理建议

### 处理完整数据集
```bash
# 建议分stage处理
for stage in E70 E80 EX05; do
  echo "处理 $stage..."
  python build_stage_h5ad_mps.py \
    --stage-dir /path/to/$stage \
    --obs-xlsx /path/to/metadata.xlsx \
    --output /path/to/${stage}_dist.h5ad \
    --batch-size 16 \
    --distance-features
done
```

### 性能优化参数

| 参数 | 推荐值 | 说明 |
|------|--------|------|
| `--batch-size` | 8-16 | MPS: 8, CUDA: 16-32 |
| `--distance-features` | 启用 | 减少99.9%特征数 |
| `--max-cells` | 测试用 | 先用50个细胞测试 |

## 🔧 环境安装

### MPS版本 (Apple Silicon)
```bash
# 确保PyTorch支持MPS
pip install torch torchvision torchaudio
python -c "import torch; print('MPS可用:', torch.backends.mps.is_available())"
```

### CUDA版本 (NVIDIA GPU)
```bash
# 根据CUDA版本安装CuPy
pip install cupy-cuda11x  # CUDA 11.x
# 或
pip install cupy-cuda12x  # CUDA 12.x
```

## 📊 输出文件说明

### 距离特征版本 (.h5ad)
- **X**: 形状 `(n_cells, 4000)`，每行一个细胞
- **obs**: 细胞metadata (cell_id, stage等)
- **var**: 特征名称 (chr1_dist1_mean, chr1_dist1_std, ...)
- **uns**: 处理参数和失败细胞列表

### 文件大小预估
- **距离特征**: ~400KB/100细胞
- **完整上三角**: ~40GB/100细胞

## 🎯 下游分析建议

### 使用scanpy分析
```python
import scanpy as sc
import anndata as ad

# 读取数据
adata = ad.read_h5ad('E70_dist.h5ad')

# 标准预处理
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 降维和聚类
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 可视化
sc.pl.umap(adata, color=['leiden', 'stage'])
```

## ⚠️ 注意事项

1. **内存监控**: 处理大数据集时监控内存使用
2. **批量大小**: 根据GPU内存调整batch_size
3. **断点续传**: 脚本支持中断恢复 (resume.json)
4. **备份数据**: 大规模处理前先备份
5. **测试先行**: 用少量细胞测试管道

## 🔗 相关文件

- `build_stage_h5ad_mps.py` - MPS加速版本 ⭐
- `build_stage_h5ad_cuda.py` - CUDA加速版本
- `build_stage_h5ad_simple.py` - CPU基础版本
- `generate_h5ad_from_cells.py` - 增量写入版本

## 📞 故障排除

### 常见问题
1. **MPS不可用**: 检查PyTorch版本和硬件支持
2. **内存不足**: 减小batch_size或使用distance_features
3. **处理中断**: 检查resume.json并重新运行
4. **速度慢**: 启用distance_features模式

---
**最后更新**: 2025-08-24
**推荐配置**: MPS + 距离特征模式
