# 🚀 快速开始指南

## 📋 使用步骤

### 1. 激活环境
```bash
# 激活schicluster环境
micromamba activate schicluster

# 验证环境
echo $CONDA_DEFAULT_ENV
# 应该显示: schicluster
```

### 2. 切换到工作目录
```bash
cd /Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/plotumap
```

### 3. 运行测试

#### 方法一：使用Python测试脚本（推荐）
```bash
python run_test.py
```
然后根据菜单选择测试选项。

#### 方法二：直接运行具体测试
```bash
# 基础功能测试
python test_upper_triangle.py

# 小规模批处理测试
python process_all_stages.py --max-cells 2 --stages E95

# 核心模块测试
python hic_upper_triangle_builder.py
```

#### 方法三：使用Shell脚本
```bash
bash run_test.sh
```

## 🎯 推荐测试流程

1. **数据可用性检查** - 确认数据存在
2. **基础功能测试** - 验证核心功能
3. **小规模测试** - 处理少量数据验证流程
4. **批量处理** - 处理所有数据

## 📊 测试选项说明

### 选项1: 基础功能测试
- 检查数据可用性
- 测试单细胞处理
- 测试小规模stage处理
- 验证特征名称生成

### 选项2: 小规模批处理测试
- 处理E95 stage（细胞数较少）
- 限制每个stage 2个细胞
- 验证完整流程

### 选项3: 核心模块测试
- 单独测试矩阵构建器
- 处理E95 stage 5个细胞

### 选项4: 数据可用性检查
- 扫描所有stage目录
- 统计可用细胞数量

### 选项5: 自定义批处理测试
- 自定义细胞数量
- 自定义处理的stages

## 🔧 故障排除

### 环境问题
```bash
# 如果环境激活失败
micromamba list environments
micromamba activate schicluster

# 检查Python路径
which python
```

### 路径问题
```bash
# 确保在正确目录
pwd
# 应该显示: /Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/plotumap

# 检查数据目录
ls -la ../../output/imputed_matrices_by_stage/
```

### 权限问题
```bash
# 确保脚本可执行
chmod +x run_test.sh
chmod +x *.py
```

## 📁 预期输出目录

成功运行后会创建以下目录和文件：
```
../../output/
├── hic_upper_triangle_h5ad/          # 主要输出目录
│   ├── E95_upper_triangle.h5ad       # 测试输出
│   └── processing_report.txt         # 处理报告
└── test_upper_triangle_h5ad/         # 测试输出目录
    └── E95_upper_triangle_test.h5ad  # 基础测试输出
```

## 📋 日志文件

- `test_upper_triangle.log` - 基础测试日志
- `process_all_stages.log` - 批处理日志
- `hic_umap_pipeline.log` - 之前系统的日志

## 💡 使用建议

1. **首次使用**: 先运行选项4检查数据，然后运行选项1进行基础测试
2. **快速验证**: 使用选项2进行小规模批处理测试
3. **生产使用**: 确认测试通过后，使用完整的批处理命令

## 🚀 完整批处理命令

测试通过后，可以使用以下命令处理所有数据：
```bash
# 处理所有stages（顺序）
python process_all_stages.py

# 并行处理（推荐）
python process_all_stages.py --parallel --max-workers 2

# 处理指定stages
python process_all_stages.py --stages E70 E80 EX05

# 限制细胞数量（用于测试）
python process_all_stages.py --max-cells 10
```