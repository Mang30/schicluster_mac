#!/bin/bash
# 代码清理脚本 - 删除重复和过时的文件

echo "🧹 开始清理 plotumap 目录..."

# 要删除的重复/过时文件列表
FILES_TO_DELETE=(
    # 重复的 h5ad 构建脚本（保留 mps, cuda, simple 三个主要版本）
    "build_stage_h5ad_fast.py"           # 被 mps 版本替代
    "build_stage_h5ad_ultra_fast.py"     # 被 mps 版本替代
    "generate_h5ad_from_cells.py"        # 增量版本，复杂且有问题
    "h5ad_generation.py"                 # 早期版本
    "hic_upper_triangle_builder.py"      # 功能重复
    
    # 重复的处理器
    "fast_batch_processor.py"            # 功能被 mps 版本覆盖
    "gpu_accelerated_processor.py"       # 功能被 cuda 版本覆盖
    "process_all_stages.py"              # 简单封装，不如直接用主脚本
    
    # 测试和调试文件
    "debug_merger.py"                    # 调试用，已完成
    "expanded_test.py"                   # 测试文件
    "fixed_test.py"                      # 测试文件
    "mini_test.py"                       # 测试文件
    "quick_validation.py"                # 验证文件
    "simple_test.py"                     # 测试文件
    "test_upper_triangle.py"             # 测试文件
    "validate_visualization.py"          # 验证文件
    "run_test.py"                        # 简单测试
    
    # 过时的工具
    "cell_annotation_matching.py"        # 早期版本，功能已集成
    "correct_stage_mapping.py"           # 一次性修正工具
    "fix_resolution.py"                  # 一次性修正工具
    "matrix_merger.py"                   # 被 robust 版本替代
    "resolution_detector.py"             # 功能已集成到主脚本
    "simple_resolution_check.py"         # 简单检查工具
    "umap_visualization.py"              # 功能已集成到 visualizer
    
    # 日志文件
    "expanded_test.log"
    "fixed_test.log"
    "hic_umap_pipeline.log"
    "test_upper_triangle.log"
    
    # 重复的文档（保留最新的核心文档）
    "EXECUTION_GUIDE.md"                 # 内容已整合到 PERFORMANCE_GUIDE
    "FINAL_SUMMARY.md"                   # 临时总结
    "hires_processing_summary.md"        # 早期总结
    "README.md"                          # 用新的 README_CLEAN.md 替代
    "stage_processing_design.md"         # 设计文档，已实现
    "TROUBLESHOOTING.md"                 # 内容已整合
    "UPPER_TRIANGLE_README.md"           # 专题文档，功能已实现
    
    # 脚本文件（功能重复）
    "run_pipeline.sh"                    # 简单封装
    "run_test.sh"                        # 测试脚本
)

# 要保留的核心文件列表（确认不会误删）
KEEP_FILES=(
    # 核心处理脚本
    "build_stage_h5ad_mps.py"            # ⭐ MPS 加速版本
    "build_stage_h5ad_cuda.py"           # CUDA 加速版本  
    "build_stage_h5ad_simple.py"         # CPU 基础版本
    "test_pipeline.py"                   # 主要的测试管道
    
    # 重要模块
    "matrix_merger_robust.py"            # 稳健的矩阵合并
    "hic_anndata_builder.py"             # AnnData 构建器
    "hic_umap_visualizer.py"             # UMAP 可视化
    
    # 配置文件
    "color_mapping.json"                 # 颜色配置
    "stage_files_mapping.csv"            # 元数据映射
    
    # 核心文档
    "CLAUDE.md"                          # 用户要求保留
    "PERFORMANCE_GUIDE.md"               # 性能指南
    "QUICK_START.md"                     # 快速开始
    "README_CLEAN.md"                    # 新的说明文档
    
    # 其他
    "FigC.png"                           # 图片文件
    "__pycache__/"                       # Python 缓存目录
)

echo "📋 将删除以下文件:"
for file in "${FILES_TO_DELETE[@]}"; do
    if [ -f "$file" ]; then
        echo "  - $file"
    fi
done

echo ""
echo "💾 将保留以下核心文件:"
for file in "${KEEP_FILES[@]}"; do
    if [ -f "$file" ] || [ -d "$file" ]; then
        echo "  ✅ $file"
    fi
done

echo ""
read -p "确认删除上述文件? (y/N): " confirm

if [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]; then
    echo "🗑️  删除文件中..."
    deleted_count=0
    
    for file in "${FILES_TO_DELETE[@]}"; do
        if [ -f "$file" ]; then
            rm "$file"
            echo "已删除: $file"
            ((deleted_count++))
        elif [ -d "$file" ]; then
            rm -rf "$file"
            echo "已删除目录: $file" 
            ((deleted_count++))
        fi
    done
    
    echo "✅ 清理完成! 删除了 $deleted_count 个文件"
    echo ""
    echo "📁 清理后的目录结构:"
    ls -la *.py *.md *.json *.csv 2>/dev/null | head -20
    
else
    echo "❌ 清理已取消"
fi
