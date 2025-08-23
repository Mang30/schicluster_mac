#!/bin/bash
# Task 6: EX15 (后期阶段独立处理)
# M4 MacBook 24GB 优化版本

echo "🚀 启动 Task 6: EX15 插补任务"
echo "💻 硬件: M4 MacBook 24GB"
echo "📊 预计细胞数: ~1,122 个"
echo "⏱️  预计时间: 3-6 小时"
echo "================================"

# 激活环境
source ~/.bashrc
micromamba activate schicluster

# 设置任务参数
TASK_NAME="Task6_EX15"
STAGES="EX15"
BATCH_SIZE=10  # M4 24GB 标准批次
LOG_FILE="task6_imputation.log"

# 创建任务日志
echo "$(date): Task 6 开始" > ../logs/$LOG_FILE

# 执行插补
echo "🔧 执行命令: python ../core/batch_mps_imputation_v2.py --stages $STAGES --batch-size $BATCH_SIZE --resolution 100000"
python ../core/batch_mps_imputation_v2.py \
    --stages $STAGES \
    --batch-size $BATCH_SIZE \
    --resolution 100000 \
    2>&1 | tee -a ../logs/$LOG_FILE

# 检查执行结果
EXIT_CODE=${PIPESTATUS[0]}
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Task 6 (EX15) 执行成功" | tee -a ../logs/$LOG_FILE
else
    echo "❌ Task 6 (EX15) 执行失败，退出码: $EXIT_CODE" | tee -a ../logs/$LOG_FILE
fi

echo "$(date): Task 6 完成，退出码: $EXIT_CODE" >> ../logs/$LOG_FILE
echo "📄 日志文件: ../logs/$LOG_FILE"

# 统计输出文件
OUTPUT_DIR="/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage"
if [ -d "$OUTPUT_DIR" ]; then
    EX15_COUNT=$(find "$OUTPUT_DIR/EX15" -name "*.npz" 2>/dev/null | wc -l)
    echo "📈 Task 6 输出统计:" | tee -a ../logs/$LOG_FILE
    echo "   EX15: $EX15_COUNT 个矩阵文件" | tee -a ../logs/$LOG_FILE
fi

exit $EXIT_CODE