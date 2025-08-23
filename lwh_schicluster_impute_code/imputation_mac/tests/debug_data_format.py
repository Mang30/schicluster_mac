#!/usr/bin/env python3
"""
测试正确的参数设置来读取您的数据格式
"""

import pandas as pd
import sys
sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')

# 测试数据读取
contact_file = "/Volumes/SumSung500/CSU/0_HiRES/converted_matrices_by_stage_20k/E70/GasdE701001/GasdE701001_chr1.txt"

print("🔍 分析数据格式...")
print("\n1. 原始数据前几行:")
with open(contact_file, 'r') as f:
    for i, line in enumerate(f):
        if i < 5:
            print(f"   {line.strip()}")

print("\n2. 用 pandas 读取:")
df = pd.read_csv(contact_file, sep='\t', header=None, index_col=None, comment='#')
print(f"   列数: {len(df.columns)}")
print(f"   前几行:")
print(df.head())

print(f"\n3. 数据范围:")
print(f"   第0列 (bin1) 范围: {df[0].min()} - {df[0].max()}")
print(f"   第1列 (bin2) 范围: {df[1].min()} - {df[1].max()}")
print(f"   第2列 (count) 范围: {df[2].min()} - {df[2].max()}")

print(f"\n4. 您的数据格式是: bin1  bin2  count")
print(f"   正确的参数应该是:")
print(f"   - 如果函数期望选择4列: [[0, 0, 1, 1]]")
print(f"   - 或者需要修改函数逻辑以处理这种格式")

# 尝试不同的参数组合
print(f"\n5. 测试不同的参数组合:")

test_params = [
    {"chrom1": 0, "pos1": 0, "chrom2": 1, "pos2": 1, "desc": "bin1, bin1, bin2, bin2"},
    {"chrom1": 0, "pos1": 1, "chrom2": 0, "pos2": 1, "desc": "bin1, bin2, bin1, bin2"},
]

for i, params in enumerate(test_params):
    try:
        print(f"\n   测试 {i+1}: {params['desc']}")
        cols = [params['chrom1'], params['pos1'], params['chrom2'], params['pos2']]
        test_df = df[cols]
        print(f"   ✅ 成功选择列: {cols}")
        print(f"   数据形状: {test_df.shape}")
    except Exception as e:
        print(f"   ❌ 失败: {e}")

print(f"\n💡 建议:")
print(f"   您的数据已经是 bin 格式，不需要染色体过滤")
print(f"   需要修改函数调用方式或使用专门的处理逻辑")
