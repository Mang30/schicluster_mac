# M4 MacBook MPS 加速单细胞 Hi-C 插补执行计划

## 硬件配置
- **设备**: M4 MacBook (Apple Silicon)
- **运行内存**: 24GB 统一内存
- **GPU加速**: MPS (Metal Performance Shaders)
- **预期性能**: 10-15x CPU 加速比

## 项目概述
对 `output/pairs2matrix_output/` 中 7 个发育阶段的单细胞 Hi-C 数据进行 MPS 加速插补处理。

## 数据统计
- **数据源**: `/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output/`
- **发育阶段**: E70, E75, E80, E85, E95, EX05, EX15 (共7个阶段)
- **总细胞数**: ~7,476个细胞
- **每细胞染色体**: 20个 (chr1-19, chrX)
- **预计总矩阵数**: ~149,520个

### 各阶段细胞统计
| 发育阶段 | 细胞数量 | 预计矩阵数 |
|---------|---------|-----------|
| E70     | ~557    | ~11,140   |
| E75     | ~1,870  | ~37,400   |
| E80     | ~559    | ~11,180   |
| E85     | ~1,125  | ~22,500   |
| E95     | ~1,108  | ~22,160   |
| EX05    | ~1,128  | ~22,560   |
| EX15    | ~1,122  | ~22,440   |
| **总计** | **7,469** | **149,380** |

## 并行任务分配策略

### 6个 Bash 任务 (M4 24GB 内存优化)
- **Task 1**: E70 + E80 (~1,116细胞) - 早期阶段合并
- **Task 2**: E75 (~1,870细胞) - 最大阶段独立处理  
- **Task 3**: E85 (~1,125细胞) - 中期阶段独立
- **Task 4**: E95 (~1,108细胞) - 中期阶段独立
- **Task 5**: EX05 (~1,128细胞) - 后期阶段独立
- **Task 6**: EX15 (~1,122细胞) - 后期阶段独立

## 技术架构

### MPS 加速配置
- **批次大小**: 8-12 (充分利用24GB内存)
- **MPS设备**: `torch.device('mps')`
- **内存管理**: 动态分配，避免内存碎片
- **温度监控**: M4芯片温度保护机制

### 输入输出路径
```
输入: /Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output/
      ├── E70/GasdE701001/GasdE701001_chr1.txt
      ├── E75/GasaE751001/GasaE751001_chr1.txt
      └── ...

输出: /Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage/
      ├── E70/GasdE701001/GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz
      ├── E75/GasaE751001/GasaE751001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz
      └── ...
```

## 执行计划

### 阶段1: 脚本开发和适配 ✅
1. ✅ 创建 CLAUDE.md 文档
2. ✅ 创建 M4 优化的 MPS 插补主控脚本
3. ✅ 修改 batch_mps_imputation.py 适配 pairs2matrix_output 数据源
4. ✅ 创建 6 个并行 bash 任务脚本
5. ✅ 配置 M4 24GB 内存优化参数
6. ✅ 更新路径从 lwh_code 到 lwh_schicluster_impute_code

### 阶段2: 监控和日志系统 ✅
1. ✅ 创建进度监控脚本
2. ✅ 设置实时日志记录
3. ✅ 系统资源监控
4. ✅ 创建环境测试脚本

### 阶段3: 测试和验证 🔄
1. 🔄 环境依赖测试
2. 🔄 小规模测试（单细胞验证）
3. 🔄 完整并行执行
4. 🔄 质量验证和报告

## 性能预估 (M4 24GB 配置)

### 处理时间预估
- **单染色体时间**: ~0.08-0.15秒 (M4 优化)
- **单细胞时间**: ~1.6-3.0秒 (20个染色体)
- **Task 1 (E70+E80)**: ~3-5小时
- **Task 2 (E75)**: ~8-12小时  
- **Task 3-6**: 各 ~3-6小时
- **总预计时间**: 12-20小时 (6任务并行)

### 资源使用
- **CPU 利用率**: 80-95% (6个Python进程)
- **GPU 利用率**: 85-95% (MPS并发)
- **内存使用**: 18-22GB (充分利用24GB)
- **温度控制**: 自动频率调节保护M4芯片

## 质量保证

### 数据验证
- **输入格式**: bin1 bin2 count (制表符分隔)
- **输出格式**: .npz 插补后的接触矩阵
- **完整性检查**: 20个染色体文件完整性
- **数值验证**: 插补后矩阵数值范围检查

### 错误处理
- **断点续传**: 自动跳过已完成的矩阵
- **错误重试**: 失败任务自动重试机制
- **日志记录**: 详细的错误和警告记录
- **资源监控**: 内存溢出和温度异常保护

## 预期输出结构
```
output/imputed_matrices_by_stage/
├── E70/
│   ├── GasdE701001/
│   │   ├── GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz
│   │   ├── GasdE701001_chr2_pad1_std1.0_rp0.5_sqrtvc.npz
│   │   └── ... (20个染色体)
│   └── ...
├── E75/
├── E80/
├── E85/
├── E95/
├── EX05/
└── EX15/
```

## 后续应用
插补完成的 .npz 矩阵文件可直接用于：
1. 单细胞 Hi-C 数据分析
2. 染色质结构研究
3. 发育阶段比较分析
4. schicluster 下游分析流程

## 使用指南

### 快速开始
```bash
# 进入脚本目录
cd /Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac

# 运行快速测试（自动激活schicluster环境）
bash quick_test.sh
```

### 手动步骤
```bash
# 1. 激活环境
eval "$(micromamba shell hook --shell bash)"
micromamba activate schicluster

# 2. 环境测试
python test_environment.py

# 3. 单细胞测试
python batch_mps_imputation_v2.py --stages E70 --max-cells 1

# 4. 完整执行
bash start_parallel_imputation.sh
```

### 监控命令
```bash
# 实时监控
python monitor_progress.py

# 一次性状态
python monitor_progress.py --once

# 生成报告  
python monitor_progress.py --summary
```

### 文件清单
- `batch_mps_imputation_v2.py`: M4优化的MPS插补脚本
- `parallel_imputation_controller.py`: 6任务并行控制器
- `task1-6_imputation.sh`: 6个并行任务脚本
- `monitor_progress.py`: 进度监控脚本
- `test_environment.py`: 环境测试脚本
- `quick_test.sh`: 快速测试脚本
- `start_parallel_imputation.sh`: 交互启动脚本

---
**创建时间**: 2025-08-23  
**硬件平台**: M4 MacBook 24GB  
**环境要求**: micromamba schicluster  
**预计完成**: 12-20小时 (6任务并行)  
**状态**: 开发完成，准备测试 ✅