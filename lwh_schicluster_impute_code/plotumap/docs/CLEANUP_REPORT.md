# 📁 PlotUMAP 代码整理完成报告

## 🎯 保留的核心文件

### ⭐ 主要处理脚本
- **build_stage_h5ad_mps.py** - MPS GPU 加速版本 (推荐)
- **build_stage_h5ad_cuda.py** - CUDA GPU 加速版本
- **build_stage_h5ad_simple.py** - CPU 基础版本
- **test_pipeline.py** - 主要的测试管道

### 🔧 核心模块
- **matrix_merger_robust.py** - 稳健的矩阵合并工具
- **hic_anndata_builder.py** - AnnData 构建器
- **hic_umap_visualizer.py** - UMAP 可视化工具

### 📋 配置文件
- **color_mapping.json** - 颜色配置
- **stage_files_mapping.csv** - 元数据映射

### 📚 文档
- **CLAUDE.md** - 用户对话记录 (用户要求保留)
- **PERFORMANCE_GUIDE.md** - 性能指南和使用说明
- **QUICK_START.md** - 快速开始指南
- **README_CLEAN.md** - 整理后的说明文档

### 🖼️ 其他
- **FigC.png** - 参考图片
- **__pycache__/** - Python 缓存目录

## 🗑️ 删除的重复/过时文件 (共计 31 个)

### 重复的构建脚本 (5个)
- build_stage_h5ad_fast.py
- build_stage_h5ad_ultra_fast.py
- generate_h5ad_from_cells.py
- h5ad_generation.py
- hic_upper_triangle_builder.py

### 重复的处理器 (3个)
- fast_batch_processor.py
- gpu_accelerated_processor.py
- process_all_stages.py

### 测试和调试文件 (10个)
- debug_merger.py
- expanded_test.py
- fixed_test.py
- mini_test.py
- quick_validation.py
- simple_test.py
- test_upper_triangle.py
- validate_visualization.py
- run_test.py
- run_test.sh

### 过时工具 (7个)
- cell_annotation_matching.py
- correct_stage_mapping.py
- fix_resolution.py
- matrix_merger.py
- resolution_detector.py
- simple_resolution_check.py
- umap_visualization.py

### 重复文档 (8个)
- EXECUTION_GUIDE.md
- FINAL_SUMMARY.md
- hires_processing_summary.md
- README.md
- stage_processing_design.md
- TROUBLESHOOTING.md
- UPPER_TRIANGLE_README.md
- run_pipeline.sh

### 日志文件 (4个)
- expanded_test.log
- fixed_test.log
- hic_umap_pipeline.log
- test_upper_triangle.log

## 📊 清理前后对比

- **清理前**: 40+ 个文件，大量重复和过时代码
- **清理后**: 13 个核心文件，功能明确，易于维护

## 🚀 推荐使用流程

1. **快速开始**: 阅读 `QUICK_START.md`
2. **性能调优**: 查看 `PERFORMANCE_GUIDE.md`
3. **主要脚本**: 使用 `build_stage_h5ad_mps.py` (Mac M1/M2) 或 `build_stage_h5ad_cuda.py` (NVIDIA GPU)
4. **可视化**: 使用 `hic_umap_visualizer.py` 生成 UMAP 图
5. **问题排查**: 参考 `CLAUDE.md` 中的对话历史

## ✅ 代码整理完成

所有重复和过时的代码已被识别，核心功能保留在精简的文件结构中。用户要求保留的 `CLAUDE.md` 文件已确认保留。
