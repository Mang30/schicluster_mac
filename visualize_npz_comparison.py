#!/usr/bin/env python3
"""
可视化对比 npz 文件内容
"""

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_sparse_matrix_from_npz(filepath):
    """从 npz 文件加载稀疏矩阵"""
    data = np.load(filepath, allow_pickle=True)
    
    format_type = data['format'].item().decode('utf-8')
    shape = tuple(data['shape'])
    matrix_data = data['data']
    
    if format_type == 'coo':
        row = data['row']
        col = data['col']
        matrix = sp.coo_matrix((matrix_data, (row, col)), shape=shape)
    elif format_type == 'csr':
        indices = data['indices']
        indptr = data['indptr']
        matrix = sp.csr_matrix((matrix_data, indices, indptr), shape=shape)
    else:
        raise ValueError(f"Unsupported format: {format_type}")
    
    data.close()
    return matrix

def create_comparison_plots():
    """创建对比图表"""
    data_dir = Path("/Volumes/SumSung500/CSU/0_HiRES/data/example_data")
    npz_files = list(data_dir.glob("*.npz"))
    
    # 加载所有矩阵
    matrices = {}
    for npz_file in npz_files:
        matrices[npz_file.stem] = load_sparse_matrix_from_npz(str(npz_file))
    
    # 创建多子图布局
    fig = plt.figure(figsize=(20, 16))
    
    # 1. 矩阵大小和稀疏度对比
    ax1 = plt.subplot(3, 3, 1)
    names = list(matrices.keys())
    shapes = [f"{m.shape[0]}×{m.shape[1]}" for m in matrices.values()]
    sparsities = [(1 - m.nnz / (m.shape[0] * m.shape[1])) * 100 for m in matrices.values()]
    
    x_pos = np.arange(len(names))
    bars = ax1.bar(x_pos, sparsities, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
    ax1.set_xlabel('文件')
    ax1.set_ylabel('稀疏度 (%)')
    ax1.set_title('矩阵稀疏度对比')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([name.replace('_', '\n') for name in names], rotation=45, ha='right')
    
    # 在每个柱子上添加形状信息
    for i, (bar, shape) in enumerate(zip(bars, shapes)):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{shape}', ha='center', va='bottom', fontsize=8)
    
    # 2. 非零元素数量对比
    ax2 = plt.subplot(3, 3, 2)
    nnz_counts = [m.nnz for m in matrices.values()]
    bars2 = ax2.bar(x_pos, nnz_counts, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
    ax2.set_xlabel('文件')
    ax2.set_ylabel('非零元素数')
    ax2.set_title('非零元素数量对比')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([name.replace('_', '\n') for name in names], rotation=45, ha='right')
    ax2.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    
    # 3. 数据类型和数值范围
    ax3 = plt.subplot(3, 3, 3)
    dtypes = [str(m.dtype) for m in matrices.values()]
    ranges = []
    for m in matrices.values():
        if m.nnz > 0:
            ranges.append(f"[{np.min(m.data):.2f}, {np.max(m.data):.2f}]")
        else:
            ranges.append("[空]")
    
    # 创建表格显示数据类型和范围
    table_data = []
    for i, name in enumerate(names):
        table_data.append([name.replace('_', '\n'), dtypes[i], ranges[i]])
    
    ax3.axis('tight')
    ax3.axis('off')
    table = ax3.table(cellText=table_data,
                     colLabels=['文件名', '数据类型', '数值范围'],
                     cellLoc='center',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.2, 1.5)
    ax3.set_title('数据类型和数值范围')
    
    # 4-7. 每个矩阵的热图（小样本）
    for i, (name, matrix) in enumerate(matrices.items()):
        ax = plt.subplot(3, 3, 4 + i)
        
        # 取样本进行可视化
        if matrix.shape[0] > 0 and matrix.shape[1] > 0:
            sample_size = min(50, matrix.shape[0], matrix.shape[1])
            if sp.isspmatrix_coo(matrix):
                sample = matrix.tocsr()[:sample_size, :sample_size].toarray()
            else:
                sample = matrix[:sample_size, :sample_size].toarray()
            
            im = ax.imshow(sample, cmap='viridis', aspect='auto')
            ax.set_title(f'{name.replace("_", " ")}\n({sample_size}×{sample_size} 样本)', fontsize=8)
            ax.set_xlabel('列索引')
            ax.set_ylabel('行索引')
            plt.colorbar(im, ax=ax, shrink=0.6)
        else:
            ax.text(0.5, 0.5, '空矩阵', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(name.replace("_", " "))
    
    # 8. 数值分布对比（直方图）
    ax8 = plt.subplot(3, 3, 8)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    for i, (name, matrix) in enumerate(matrices.items()):
        if matrix.nnz > 0:
            # 对于大数据集，随机采样
            data = matrix.data
            if len(data) > 10000:
                sample_indices = np.random.choice(len(data), 10000, replace=False)
                data = data[sample_indices]
            
            ax8.hist(data, bins=50, alpha=0.6, label=name.replace('_', ' '), 
                    color=colors[i % len(colors)], density=True)
    
    ax8.set_xlabel('数值')
    ax8.set_ylabel('密度')
    ax8.set_title('数值分布对比')
    ax8.legend(fontsize=8)
    ax8.set_yscale('log')
    
    # 9. 文件大小对比
    ax9 = plt.subplot(3, 3, 9)
    file_sizes = []
    for npz_file in npz_files:
        size_mb = npz_file.stat().st_size / (1024 * 1024)
        file_sizes.append(size_mb)
    
    bars3 = ax9.bar(x_pos, file_sizes, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
    ax9.set_xlabel('文件')
    ax9.set_ylabel('文件大小 (MB)')
    ax9.set_title('文件大小对比')
    ax9.set_xticks(x_pos)
    ax9.set_xticklabels([name.replace('_', '\n') for name in names], rotation=45, ha='right')
    
    # 在柱子上添加具体数值
    for bar, size in zip(bars3, file_sizes):
        height = bar.get_height()
        ax9.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{size:.1f}MB', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig('/Volumes/SumSung500/CSU/0_HiRES/npz_files_comparison.png', 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    print("可视化图表已保存为: npz_files_comparison.png")

if __name__ == "__main__":
    create_comparison_plots()
