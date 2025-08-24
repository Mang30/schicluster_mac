# Hi-C 单细胞数据处理工具集

## 📁 核心处理脚本

### 主要工具 (推荐使用)
- `build_stage_h5ad_mps.py` - **MPS GPU加速版本** ⭐ (Apple Silicon)
- `build_stage_h5ad_cuda.py` - **CUDA GPU加速版本** (NVIDIA GPU)
- `build_stage_h5ad_simple.py` - **CPU基础版本** (备用)

### 辅助模块
- `matrix_merger_robust.py` - 稳健的矩阵合并工具
- `hic_anndata_builder.py` - AnnData构建器
- `hic_umap_visualizer.py` - UMAP可视化工具

## 📋 配置文件
- `color_mapping.json` - 可视化颜色配置
- `stage_files_mapping.csv` - 样本元数据映射

## 📖 文档
- `PERFORMANCE_GUIDE.md` - **性能对比与使用指南** ⭐
- `README.md` - 项目说明
- `QUICK_START.md` - 快速开始指南

## 🚀 快速使用

### MPS加速版本 (推荐)
```bash
python build_stage_h5ad_mps.py \
  --stage-dir /path/to/stage \
  --obs-xlsx /path/to/metadata.xlsx \
  --output /path/to/output.h5ad \
  --distance-features
```

### 处理所有阶段
```bash
for stage in E70 E80 EX05; do
  python build_stage_h5ad_mps.py \
    --stage-dir /path/to/$stage \
    --obs-xlsx /path/to/metadata.xlsx \
    --output /path/to/${stage}.h5ad \
    --distance-features
done
```

## 📊 性能对比

| 方法 | 特征维度 | 3细胞耗时 | 推荐度 |
|------|----------|-----------|--------|
| MPS距离特征 | 4000维 | 1.5分钟 | ⭐⭐⭐ |
| CUDA完整 | 1.15亿维 | 2-3分钟 | ⭐⭐ |
| CPU基础 | 1.15亿维 | >5分钟 | ⭐ |

---
**维护者**: AI Assistant  
**最后更新**: 2025-08-24
