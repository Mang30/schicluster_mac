#!/bin/bash
# Task 5: EX05 (后期阶段独立处理)
# M4 MacBook 24GB 优化版本

echo "🚀 启动 Task 5: EX05 插补任务"
echo "💻 硬件: M4 MacBook 24GB"
echo "📊 预计细胞数: ~1,128 个"
echo "⏱️  预计时间: 3-6 小时"
echo "================================"

# 激活环境
source ~/.bashrc
micromamba activate schicluster

# 设置任务参数
TASK_NAME="Task5_EX05"
STAGES="EX05"
BATCH_SIZE=10  # M4 24GB 标准批次
LOG_FILE="task5_imputation.log"

# 创建任务日志
echo "$(date): Task 5 开始" > ../logs/$LOG_FILE

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
    echo "✅ Task 5 (EX05) 执行成功" | tee -a ../logs/$LOG_FILE
else
    echo "❌ Task 5 (EX05) 执行失败，退出码: $EXIT_CODE" | tee -a ../logs/$LOG_FILE
fi

echo "$(date): Task 5 完成，退出码: $EXIT_CODE" >> ../logs/$LOG_FILE
echo "📄 日志文件: ../logs/$LOG_FILE"

# 统计输出文件
OUTPUT_DIR="/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage"
if [ -d "$OUTPUT_DIR" ]; then
    EX05_COUNT=$(find "$OUTPUT_DIR/EX05" -name "*.npz" 2>/dev/null | wc -l)
    echo "📈 Task 5 输出统计:" | tee -a ../logs/$LOG_FILE
    echo "   EX05: $EX05_COUNT 个矩阵文件" | tee -a ../logs/$LOG_FILE
fi

exit $EXIT_CODE