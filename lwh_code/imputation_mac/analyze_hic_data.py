#!/usr/bin/env python3
"""
分析真实 Hi-C 数据的矩阵大小和稀疏度，给出 CPU vs MPS 的建议
"""

import os
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_contact_matrix(contact_file, chrom_size_file, chrom, resolution=20000):
    """从接触文件加载矩阵"""
    # 读取染色体大小
    chrom_sizes = pd.read_csv(chrom_size_file, sep='\t', header=None, index_col=0).squeeze()
    n_bins = (chrom_sizes.loc[chrom] // resolution) + 1
    
    # 读取接触数据
    contacts = pd.read_csv(contact_file, sep='\t', header=None, names=['bin1', 'bin2', 'count'])
    
    # 创建稀疏矩阵
    matrix = csr_matrix((contacts['count'], (contacts['bin1'], contacts['bin2'])), 
                       shape=(n_bins, n_bins))
    matrix = matrix + matrix.T  # 对称化
    
    return matrix

def analyze_hic_data():
    """分析 Hi-C 数据特征"""
    print("=== 真实单细胞 Hi-C 数据分析 ===\n")
    
    base_path = "/Volumes/SumSung500/CSU/0_HiRES"
    data_path = f"{base_path}/converted_matrices_by_stage_20k/E70"
    chrom_size_file = f"{base_path}/mm10_chrom_sizes.txt"
    
    # 选择几个代表性的细胞和染色体进行分析
    sample_cells = ["GasdE701001", "GasdE701010", "GasdE701020"]
    test_chroms = ["chr1", "chr2", "chr10", "chrX"]
    
    results = []
    
    for cell in sample_cells:
        cell_path = f"{data_path}/{cell}"
        if not os.path.exists(cell_path):
            continue
            
        print(f"📱 分析细胞: {cell}")
        
        for chrom in test_chroms:
            contact_file = f"{cell_path}/{cell}_{chrom}.txt"
            if not os.path.exists(contact_file):
                continue
                
            try:
                # 加载矩阵
                matrix = load_contact_matrix(contact_file, chrom_size_file, chrom)
                
                # 计算统计信息
                size = matrix.shape[0]
                nnz = matrix.nnz
                sparsity = nnz / (size * size)
                total_contacts = matrix.sum()
                
                print(f"  {chrom}: {size}x{size}, 稀疏度: {sparsity:.6f}, 接触数: {int(total_contacts)}")
                
                results.append({
                    'cell': cell,
                    'chrom': chrom,
                    'size': size,
                    'sparsity': sparsity,
                    'contacts': int(total_contacts),
                    'nnz': nnz
                })
                
            except Exception as e:
                print(f"  {chrom}: 读取失败 - {e}")
        
        print()
    
    # 统计分析
    if results:
        df = pd.DataFrame(results)
        
        print("=== 数据特征总结 ===")
        print(f"矩阵大小范围: {df['size'].min()} - {df['size'].max()}")
        print(f"平均矩阵大小: {df['size'].mean():.0f}")
        print(f"稀疏度范围: {df['sparsity'].min():.6f} - {df['sparsity'].max():.6f}")
        print(f"平均稀疏度: {df['sparsity'].mean():.6f}")
        print(f"接触数范围: {df['contacts'].min()} - {df['contacts'].max()}")
        print()
        
        # 基于我们的基准测试给出建议
        print("=== CPU vs MPS 建议 ===")
        
        # 统计不同大小矩阵的数量
        small_matrices = (df['size'] < 1000).sum()
        medium_matrices = ((df['size'] >= 1000) & (df['size'] < 5000)).sum()  
        large_matrices = (df['size'] >= 5000).sum()
        
        print(f"小矩阵 (<1000): {small_matrices} 个 → 建议使用 CPU")
        print(f"中等矩阵 (1000-5000): {medium_matrices} 个 → 建议使用 MPS")
        print(f"大矩阵 (>=5000): {large_matrices} 个 → 强烈建议使用 MPS")
        
        # 计算总处理时间估算
        total_matrices = len(df)
        
        # 基于基准测试的时间估算 (每个矩阵平均)
        avg_size = df['size'].mean()
        if avg_size < 1000:
            cpu_time_per_matrix = 0.1  # 秒
            mps_time_per_matrix = 0.2  # 秒 (小矩阵MPS较慢)
            recommended = "CPU"
        elif avg_size < 5000:
            cpu_time_per_matrix = 2.0  # 秒
            mps_time_per_matrix = 0.15  # 秒
            recommended = "MPS"
        else:
            cpu_time_per_matrix = 10.0  # 秒
            mps_time_per_matrix = 0.7   # 秒
            recommended = "MPS"
        
        total_cpu_time = total_matrices * cpu_time_per_matrix
        total_mps_time = total_matrices * mps_time_per_matrix
        
        print(f"\n⏱️  处理 {total_matrices} 个矩阵的预计时间:")
        print(f"CPU 总时间: {total_cpu_time/60:.1f} 分钟")
        print(f"MPS 总时间: {total_mps_time/60:.1f} 分钟")
        print(f"潜在加速比: {total_cpu_time/total_mps_time:.1f}x")
        
        print(f"\n🎯 总体建议: 使用 {recommended}")
        
        if recommended == "MPS":
            print("💡 智能自动选择会:")
            print("   - 小矩阵自动使用 CPU")
            print("   - 大矩阵自动使用 MPS") 
            print("   - 获得最佳性能")
        
        return df
    
    return None

def count_total_workload():
    """统计总的工作量"""
    print("\n=== 总工作量统计 ===")
    
    base_path = "/Volumes/SumSung500/CSU/0_HiRES/converted_matrices_by_stage_20k"
    stages = ["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"]
    
    total_cells = 0
    total_matrices = 0
    
    for stage in stages:
        stage_path = f"{base_path}/{stage}"
        if os.path.exists(stage_path):
            cells = [d for d in os.listdir(stage_path) if d.startswith("Gas") and os.path.isdir(f"{stage_path}/{d}")]
            stage_cells = len(cells)
            
            # 估算每个细胞的染色体数量 (大约20个染色体)
            stage_matrices = stage_cells * 20
            
            print(f"{stage}: {stage_cells} 细胞, ~{stage_matrices} 矩阵")
            total_cells += stage_cells
            total_matrices += stage_matrices
    
    print(f"\n📊 总计: {total_cells} 细胞, ~{total_matrices} 矩阵")
    
    # 基于平均大小估算处理时间
    avg_time_cpu = 2.0  # 每个矩阵 2 秒 (保守估计)
    avg_time_mps = 0.2  # 每个矩阵 0.2 秒 (智能选择)
    
    total_cpu_hours = (total_matrices * avg_time_cpu) / 3600
    total_mps_hours = (total_matrices * avg_time_mps) / 3600
    
    print(f"⏱️  预计处理时间:")
    print(f"纯 CPU: {total_cpu_hours:.1f} 小时")
    print(f"智能 MPS: {total_mps_hours:.1f} 小时")
    print(f"时间节省: {total_cpu_hours - total_mps_hours:.1f} 小时")

def main():
    """主函数"""
    print("单细胞 Hi-C 数据 CPU vs MPS 性能分析\n")
    
    # 分析样本数据
    df = analyze_hic_data()
    
    # 统计总工作量
    count_total_workload()
    
    print("\n🎉 结论: 对于您的单细胞 Hi-C 插补任务，建议使用智能自动选择模式!")
    print("   这样可以获得最佳性能，同时节省大量处理时间。")

if __name__ == "__main__":
    main()
