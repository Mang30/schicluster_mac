#!/usr/bin/env python3
"""
MPS 加速 Hi-C 插补使用指南和快速开始脚本
"""

print("🚀 MPS 加速单细胞 Hi-C 插补 - 使用指南")
print("=" * 60)

print("\n📋 使用方法:")
print("1. 测试单个 stage (推荐先试运行):")
print("   micromamba activate schicluster")
print("   python batch_mps_imputation.py --stages E70 --max-cells 2")
print()

print("2. 处理特定 stages:")
print("   python batch_mps_imputation.py --stages E70 E75")
print()

print("3. 处理所有主要 stages:")
print("   python batch_mps_imputation.py --stages E70 E75 E80")
print()

print("4. 如果需要使用 CPU (调试):")
print("   python batch_mps_imputation.py --stages E70 --disable-mps")
print()

print("📁 输出结构:")
print("   /Volumes/SumSung500/CSU/0_HiRES/imputed_matrices_by_stage_20k/")
print("   ├── E70/")
print("   │   ├── GasdE701001/")
print("   │   │   ├── GasdE701001_chr1_pad1_std1.0_rp0.5_sqrtvc.npz")
print("   │   │   ├── GasdE701001_chr2_pad1_std1.0_rp0.5_sqrtvc.npz")
print("   │   │   └── ...")
print("   │   └── ...")
print("   ├── E75/")
print("   └── ...")
print()

print("📊 处理统计:")
print("   - E70: ~557 细胞 x 20 染色体 = ~11,140 矩阵")
print("   - E75: ~1,870 细胞 x 20 染色体 = ~37,400 矩阵") 
print("   - E80: ~559 细胞 x 20 染色体 = ~11,180 矩阵")
print("   - 总计: ~59,720 矩阵")
print()

print("⏱️  预计时间 (MPS 加速):")
print("   - E70: ~1.8 小时")
print("   - E75: ~6.2 小时") 
print("   - E80: ~1.9 小时")
print("   - 总计: ~10 小时 (vs 纯 CPU 的 ~100 小时)")
print()

print("🔍 监控和日志:")
print("   - 实时日志输出到终端")
print("   - 详细日志保存到 'mps_imputation.log'")
print("   - 每个细胞完成后显示进度和预计剩余时间")
print()

print("💡 推荐的处理策略:")
print("   1. 先运行小规模测试确保一切正常")
print("   2. 从较小的 stage 开始 (E70 或 E80)")
print("   3. 如果有问题可以随时 Ctrl+C 停止")
print("   4. 支持断点续传 - 已完成的矩阵会自动跳过")
print()

# 交互式开始
import sys

print("🎯 现在开始吗？")
print("1. 测试运行 (E70, 2个细胞)")
print("2. 处理 E70 完整")
print("3. 处理 E70 + E75")  
print("4. 处理所有主要 stages (E70 + E75 + E80)")
print("5. 查看当前进度")
print("6. 退出")

choice = input("\n请输入选择 (1-6): ").strip()

commands = {
    "1": "python batch_mps_imputation.py --stages E70 --max-cells 2",
    "2": "python batch_mps_imputation.py --stages E70", 
    "3": "python batch_mps_imputation.py --stages E70 E75",
    "4": "python batch_mps_imputation.py --stages E70 E75 E80",
    "5": "find /Volumes/SumSung500/CSU/0_HiRES/imputed_matrices_by_stage_20k -name '*.npz' | wc -l && echo '个矩阵已完成'"
}

if choice in commands:
    print(f"\n🚀 执行命令:")
    print(f"micromamba activate schicluster && {commands[choice]}")
    print(f"\n复制上面的命令到终端运行即可！")
elif choice == "6":
    print("👋 再见！")
else:
    print("❌ 无效选择，请重新运行脚本")
