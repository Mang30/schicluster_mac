# 🎯 HiRES Hi-C 数据处理工具包

## 📁 文件结构

```
plotumap/
├── 🚀 核心处理脚本
│   ├── build_stage_h5ad_mps.py      # MPS GPU 加速版本 (推荐 Apple Silicon)
│   ├── build_stage_h5ad_cuda.py     # CUDA GPU 加速版本 (NVIDIA GPU)
│   ├── build_stage_h5ad_simple.py   # CPU 基础版本
│   └── test_pipeline.py             # 测试管道
│
├── 🔧 工具模块  
│   ├── matrix_merger_robust.py      # 稳健的矩阵合并
│   ├── hic_anndata_builder.py       # AnnData 构建器
│   └── hic_umap_visualizer.py       # UMAP 可视化
│
├── ⚙️ 配置文件
│   ├── color_mapping.json           # 发育阶段颜色配置
│   └── stage_files_mapping.csv      # 细胞元数据映射
│
├── 📚 文档
│   ├── README_CLEAN.md              # 项目说明 (主要)
│   ├── PERFORMANCE_GUIDE.md         # 性能指南
│   ├── QUICK_START.md               # 快速开始
│   ├── CLAUDE.md                    # 开发对话记录
│   └── CLEANUP_REPORT.md            # 代码清理报告
│
├── 🖼️ 资源
│   └── FigC.png                     # 参考图片
│
└── 🛠️ 辅助
    └── cleanup_code.sh              # 代码清理脚本
```

## 🌟 主要功能

1. **Hi-C 数据处理**: 从单细胞 Hi-C 矩阵构建 h5ad 格式的数据文件
2. **GPU 加速**: 支持 Apple Silicon MPS 和 NVIDIA CUDA 加速
3. **特征提取**: 上三角矩阵特征提取，支持距离分箱优化
4. **可视化**: UMAP 降维和发育阶段可视化
5. **批量处理**: 支持多个发育阶段的批量处理

## 🚀 快速使用

```bash
# Apple Silicon Mac (推荐)
python build_stage_h5ad_mps.py --stage E70 --max_cells 10

# NVIDIA GPU
python build_stage_h5ad_cuda.py --stage E70 --max_cells 10

# CPU 版本
python build_stage_h5ad_simple.py --stage E70 --max_cells 10
```

## 📊 性能表现

- **CPU 版本**: 3 细胞 ~5 分钟 (完整特征)
- **MPS 版本**: 3 细胞 ~1.5 分钟 (距离特征)
- **内存优化**: 从 114M 维降至 4K 维特征

## 🎯 代码整理成果

✅ **删除了 37 个重复/过时文件**
- 5 个重复的构建脚本
- 3 个重复的处理器
- 10 个测试和调试文件
- 7 个过时工具
- 8 个重复文档
- 4 个日志文件

✅ **保留了 13 个核心文件**
- 功能明确，易于维护
- 支持多种硬件配置
- 完整的文档体系

## 💡 按需求保留的文件

根据用户要求，特别保留了 `CLAUDE.md` 文件，记录了完整的开发对话历史和技术细节。

---

*代码已整理完成，可直接投入使用！* 🎉
