# 🚀 一键启动脚本构建计划

## 📋 需求分析

### 核心功能
1. **自动发现 stages**: 扫描 `/output/imputed_matrices_by_stage` 中的所有 stage 目录
2. **完整流水线**: merge → build h5ad → generate UMAP
3. **快速测试模式**: 支持小样本测试验证功能
4. **智能选择**: 根据硬件自动选择最优处理方式（MPS/CUDA/CPU）

### 脚本参数设计
```bash
./run.sh [选项]

选项:
  --mode MODE           运行模式: all|test|stage
  --stages STAGES       指定处理的stages (逗号分隔，如 E70,E80)
  --max-cells N         每个stage最大处理细胞数 (默认: all)
  --test-cells N        测试模式下的细胞数 (默认: 3)
  --use-distance        使用距离特征 (更快，推荐)
  --skip-merge          跳过矩阵合并步骤
  --skip-h5ad           跳过h5ad构建步骤
  --skip-umap           跳过UMAP可视化步骤
  --gpu-type TYPE       GPU类型: auto|mps|cuda|cpu (默认: auto)
  --help               显示帮助信息
```

### 运行模式
1. **全量模式** (`--mode all`): 处理所有发现的stages
2. **测试模式** (`--mode test`): 快速功能验证 (3个细胞 × 3个stages)
3. **指定模式** (`--mode stage --stages E70,E80`): 处理指定stages

### 脚本结构
```
run.sh
├── 参数解析和验证
├── 环境检测 (硬件、依赖)
├── Stage 自动发现
├── 流水线执行
│   ├── 步骤1: 矩阵合并
│   ├── 步骤2: h5ad 构建
│   └── 步骤3: UMAP 可视化
├── 进度监控和日志
└── 结果总结
```

### 智能特性
1. **硬件检测**: 自动检测 Apple Silicon MPS / NVIDIA CUDA
2. **断点续传**: 检查已完成的文件，避免重复处理
3. **并行处理**: 支持多 stage 并行处理
4. **错误处理**: 优雅处理单个stage失败，继续其他stage
5. **资源监控**: 显示处理进度、时间估算、内存使用

### 输出组织
```
output/
├── merged_matrices/          # 合并后的矩阵
├── hic_h5ad_files/          # h5ad 文件
├── umap_plots/              # UMAP 图表
├── logs/                    # 运行日志
│   ├── run_YYYYMMDD_HHMMSS.log
│   └── summary.json
└── reports/                 # 处理报告
    └── processing_report.html
```

## 🔧 技术实现要点

### 1. Stage 自动发现
```bash
# 扫描所有可用的 stage 目录
find /output/imputed_matrices_by_stage -maxdepth 1 -type d -name "E*"
```

### 2. 硬件检测逻辑
```bash
detect_gpu() {
    if [[ $(uname -m) == "arm64" ]] && python -c "import torch; print(torch.backends.mps.is_available())" 2>/dev/null | grep -q "True"; then
        echo "mps"
    elif nvidia-smi >/dev/null 2>&1; then
        echo "cuda" 
    else
        echo "cpu"
    fi
}
```

### 3. 智能处理选择
- **MPS**: `build_stage_h5ad_mps.py --distance-features`
- **CUDA**: `build_stage_h5ad_cuda.py --distance-features`  
- **CPU**: `build_stage_h5ad_simple.py`

### 4. 进度监控
```bash
# 实时显示处理进度
echo "🔄 处理 Stage E70 (1/3): 合并矩阵..."
echo "⏱️  预计剩余时间: 5分钟"
echo "💾 内存使用: 2.1GB / 16GB"
```

## 📊 预期效果

### 快速测试 (--mode test)
```bash
./run.sh --mode test --test-cells 3
```
- 3个 stages × 3个细胞 = 9个样本
- 预计运行时间: 2-5分钟
- 验证完整流水线功能

### 全量处理 (--mode all)
```bash
./run.sh --mode all --use-distance
```
- 处理所有可用stages
- 使用距离特征优化性能
- 自动硬件加速

### 定制处理
```bash
./run.sh --mode stage --stages E70,E80 --max-cells 50 --gpu-type mps
```
- 仅处理 E70 和 E80
- 每个stage最多50个细胞
- 强制使用 MPS 加速

## 🎯 优化建议

1. **提示词优化**:
   - 清晰的进度指示和时间预估
   - 彩色输出增强可读性
   - 详细的错误信息和解决建议

2. **性能优化**:
   - 默认使用距离特征 (4K vs 114M 维)
   - 智能批处理大小调整
   - 内存使用监控和预警

3. **用户体验**:
   - 交互式模式选择
   - 结果可视化预览
   - 一键式环境检查

---

**下一步**: 基于此计划实现 `run.sh` 脚本 🚀
