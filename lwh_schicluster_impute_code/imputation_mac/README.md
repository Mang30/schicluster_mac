# M4 MacBook MPS 加速单细胞 Hi-C 插补系统

## 目录结构

```
imputation_mac/
├── README.md                    # 项目说明文档
├── core/                        # 核心功能模块
│   ├── batch_mps_imputation_v2.py      # M4优化的MPS插补主程序
│   ├── parallel_imputation_controller.py # 6任务并行控制器
│   ├── monitor_progress.py              # 进度监控工具
│   └── impute_bin_format.py             # 插补核心算法
├── scripts/                     # 执行脚本
│   ├── task1_imputation.sh              # Task 1: E70+E80
│   ├── task2_imputation.sh              # Task 2: E75
│   ├── task3_imputation.sh              # Task 3: E85
│   ├── task4_imputation.sh              # Task 4: E95
│   ├── task5_imputation.sh              # Task 5: EX05
│   ├── task6_imputation.sh              # Task 6: EX15
│   ├── start_parallel_imputation.sh     # 交互启动脚本
│   └── quick_test.sh                    # 快速测试脚本
├── tests/                       # 测试和验证
│   ├── test_environment.py             # 环境测试
│   ├── test_mps_acceleration.py        # MPS加速测试
│   ├── test_mps_large_matrix.py        # 大矩阵测试
│   ├── test_parallel_chroms.py         # 并行染色体测试
│   ├── benchmark_mps_vs_cpu.py         # 性能基准测试
│   ├── debug_data_format.py            # 数据格式调试
│   └── analyze_hic_data.py             # Hi-C数据分析
├── logs/                        # 日志文件
│   └── mps_imputation.log              # MPS插补日志
├── data/                        # 数据文件
│   ├── mm10_chrom_sizes.txt            # 小鼠染色体长度
│   └── test_output.npz                 # 测试输出文件
└── docs/                        # 文档
    ├── CLAUDE.md                       # 项目执行计划
    ├── README_MPS_ACCELERATION.md      # MPS加速说明
    ├── README_MPS_USAGE.md             # MPS使用指南
    └── README_PARALLEL_CHROMS.md       # 并行染色体处理说明
```

## 快速开始

### 方式一：一键启动（推荐）
```bash
# 进入项目目录
cd /Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac

# 一键启动主菜单
bash run.sh
```

### 方式二：手动步骤
```bash
# 1. 激活 schicluster 环境
eval "$(micromamba shell hook --shell bash)"
micromamba activate schicluster

# 2. 进入项目目录
cd /Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac

# 3. 选择操作
bash scripts/quick_test.sh              # 快速测试
bash scripts/start_parallel_imputation.sh  # 完整执行
python core/monitor_progress.py         # 进度监控
```

## 核心功能

### 主要程序
- **batch_mps_imputation_v2.py**: M4 MacBook 24GB 优化的MPS加速插补程序
- **parallel_imputation_controller.py**: 管理6个并行任务的主控制器
- **monitor_progress.py**: 实时监控插补进度和系统资源

### 并行策略
- **Task 1**: E70 + E80 (~1,116细胞) 
- **Task 2**: E75 (~1,870细胞) - 最大阶段
- **Task 3**: E85 (~1,125细胞)
- **Task 4**: E95 (~1,108细胞)  
- **Task 5**: EX05 (~1,128细胞)
- **Task 6**: EX15 (~1,122细胞)

## 性能特点

- **硬件优化**: 专为M4 MacBook 24GB设计
- **MPS加速**: 10-15x CPU加速比
- **并行处理**: 6任务同时执行
- **内存管理**: 智能批次调整
- **温度保护**: 自动频率调节

## 监控工具

```bash
# 实时监控
python core/monitor_progress.py

# 一次性状态查看
python core/monitor_progress.py --once

# 生成汇总报告
python core/monitor_progress.py --summary
```

## 测试工具

```bash
# 环境测试
python tests/test_environment.py

# MPS性能测试
python tests/benchmark_mps_vs_cpu.py

# 数据格式调试
python tests/debug_data_format.py
```

## 输出数据

```
/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage/
├── E70/
├── E75/
├── E80/
├── E85/
├── E95/
├── EX05/
└── EX15/
```

每个细胞生成20个染色体的.npz插补矩阵文件。

---
**硬件**: M4 MacBook 24GB  
**环境**: micromamba schicluster  
**预计时间**: 12-20小时 (6任务并行)  
**状态**: 开发完成，准备使用 ✅