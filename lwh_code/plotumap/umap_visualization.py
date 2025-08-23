#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UMAP可视化分析脚本
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import json
from collections import Counter


def load_h5ad_files(input_dir):
    """
    加载所有Stage的h5ad文件
    
    Args:
        input_dir (str): 输入目录路径
        
    Returns:
        dict: 包含所有Stage数据的字典
    """
    stage_files = [f for f in os.listdir(input_dir) if f.startswith('stage_') and f.endswith('.h5ad')]
    stage_data = {}
    
    for file in stage_files:
        stage_name = file.replace('stage_', '').replace('.h5ad', '')
        file_path = os.path.join(input_dir, file)
        print(f"加载Stage {stage_name} 数据...")
        adata = sc.read_h5ad(file_path)
        stage_data[stage_name] = adata
        print(f"Stage {stage_name} 数据加载完成，形状: {adata.shape}")
    
    return stage_data


def create_unified_color_mapping(stage_data, color_mapping_file=None):
    """
    创建统一的颜色映射方案
    
    Args:
        stage_data (dict): 所有Stage数据
        color_mapping_file (str): 颜色映射文件路径
        
    Returns:
        dict: 统一的颜色映射
    """
    # 获取所有唯一的细胞类型
    all_celltypes = set()
    for adata in stage_data.values():
        all_celltypes.update(adata.obs['Celltype'].unique())
    
    all_celltypes = sorted(list(all_celltypes))
    print(f"总共发现 {len(all_celltypes)} 种细胞类型:")
    for ct in all_celltypes:
        print(f"  - {ct}")
    
    # 创建统一的颜色映射
    if color_mapping_file and os.path.exists(color_mapping_file):
        # 从文件加载颜色映射
        with open(color_mapping_file, 'r') as f:
            color_mapping = json.load(f)
    else:
        # 生成新的颜色映射
        colors = sns.color_palette("husl", len(all_celltypes))
        color_mapping = dict(zip(all_celltypes, colors))
        
        # 保存颜色映射到文件
        if color_mapping_file:
            with open(color_mapping_file, 'w') as f:
                # 将numpy数组转换为列表以便JSON序列化
                color_mapping_serializable = {k: [float(c) for c in v] for k, v in color_mapping.items()}
                json.dump(color_mapping_serializable, f, indent=2)
            print(f"颜色映射已保存到: {color_mapping_file}")
    
    return color_mapping


def plot_umap_single_stage(adata, stage_name, output_dir, color_mapping):
    """
    绘制单个Stage的UMAP图
    
    Args:
        adata (AnnData): 包含数据的AnnData对象
        stage_name (str): Stage名称
        output_dir (str): 输出目录路径
        color_mapping (dict): 颜色映射
    """
    print(f"处理Stage {stage_name} 的UMAP可视化...")
    
    # 检查是否有足够的细胞进行UMAP分析
    if adata.n_obs < 10:
        print(f"Stage {stage_name} 的细胞数量不足(少于10个)，跳过UMAP分析")
        return
    
    # 数据预处理
    # 确保数据是密集矩阵（如果数据量不大）
    if adata.n_obs < 1000:  # 对于小数据集使用密集矩阵
        adata.X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    # 标准化数据
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 计算高度可变基因（如果适用）
    try:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        # sc.pl.highly_variable_genes(adata)
        adata = adata[:, adata.var.highly_variable]
    except:
        print(f"Stage {stage_name} 数据可能不是基因表达数据，跳过高度可变基因计算")
    
    # 运行PCA
    try:
        sc.tl.pca(adata, svd_solver='arpack')
    except:
        print(f"Stage {stage_name} PCA计算失败，使用原始数据进行UMAP")
    
    # 运行UMAP
    try:
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
    except Exception as e:
        print(f"Stage {stage_name} UMAP计算失败: {e}")
        return
    
    # 创建UMAP图
    plt.figure(figsize=(10, 8))
    
    # 获取细胞类型和对应颜色
    celltypes = adata.obs['Celltype']
    colors = [color_mapping.get(ct, (0.5, 0.5, 0.5)) for ct in celltypes]  # 默认灰色
    
    # 绘制UMAP散点图
    scatter = plt.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], 
                         c=colors, alpha=0.7, s=20)
    
    # 添加图例
    # 获取唯一的细胞类型及其颜色
    unique_celltypes = adata.obs['Celltype'].unique()
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=color_mapping.get(ct, (0.5, 0.5, 0.5)), markersize=8, 
                                 label=ct) for ct in unique_celltypes]
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.title(f'UMAP Visualization - Stage {stage_name}')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.tight_layout()
    
    # 保存图像
    output_file = os.path.join(output_dir, f"umap_{stage_name}.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Stage {stage_name} 的UMAP图已保存到: {output_file}")


def plot_umap_all_stages(stage_data, output_dir, color_mapping):
    """
    为所有Stage绘制UMAP图
    
    Args:
        stage_data (dict): 所有Stage数据
        output_dir (str): 输出目录路径
        color_mapping (dict): 颜色映射
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 按Stage分组处理
    for stage_name, adata in stage_data.items():
        print(f"处理Stage: {stage_name}")
        
        # 绘制UMAP图
        plot_umap_single_stage(adata, stage_name, output_dir, color_mapping)


def plot_umap_composite(stage_data, output_dir, color_mapping):
    """
    绘制复合UMAP图，按照FigC.png样式
    
    Args:
        stage_data (dict): 所有Stage数据
        output_dir (str): 输出目录路径
        color_mapping (dict): 颜色映射
    """
    print("生成复合UMAP图...")
    
    # 定义Stage顺序
    stage_order = ['E70', 'E75', 'E80', 'E85', 'E95', 'EX05', 'EX15']
    
    # 过滤掉不存在的Stage
    available_stages = [stage for stage in stage_order if stage in stage_data]
    
    # 创建复合图
    n_stages = len(available_stages)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    # 隐藏最后一个子图（因为我们只有7个Stage）
    axes[-1].set_visible(False)
    
    # 为每个Stage绘制UMAP图
    for i, stage_name in enumerate(available_stages):
        ax = axes[i]
        adata = stage_data[stage_name]
        
        # 检查是否有足够的细胞进行UMAP分析
        if adata.n_obs < 10:
            ax.text(0.5, 0.5, f'{stage_name}\n(细胞数不足)', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Stage {stage_name}')
            continue
        
        # 数据预处理
        if adata.n_obs < 1000:  # 对于小数据集使用密集矩阵
            adata.X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
        
        # 标准化数据
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 计算高度可变基因（如果适用）
        try:
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            adata = adata[:, adata.var.highly_variable]
        except:
            print(f"Stage {stage_name} 数据可能不是基因表达数据，跳过高度可变基因计算")
        
        # 运行PCA
        try:
            sc.tl.pca(adata, svd_solver='arpack')
        except:
            print(f"Stage {stage_name} PCA计算失败，使用原始数据进行UMAP")
        
        # 运行UMAP
        try:
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
            sc.tl.umap(adata)
        except Exception as e:
            print(f"Stage {stage_name} UMAP计算失败: {e}")
            ax.text(0.5, 0.5, f'{stage_name}\n(UMAP失败)', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Stage {stage_name}')
            continue
        
        # 获取细胞类型和对应颜色
        celltypes = adata.obs['Celltype']
        colors = [color_mapping.get(ct, (0.5, 0.5, 0.5)) for ct in celltypes]  # 默认灰色
        
        # 绘制UMAP散点图
        ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], 
                  c=colors, alpha=0.7, s=10)
        
        ax.set_title(f'Stage {stage_name}')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
    
    # 添加统一图例
    unique_celltypes = []
    legend_colors = []
    for stage_name in available_stages:
        adata = stage_data[stage_name]
        if adata.n_obs >= 10:  # 只考虑成功绘制的Stage
            uct = adata.obs['Celltype'].unique()
            for ct in uct:
                if ct not in unique_celltypes:
                    unique_celltypes.append(ct)
                    legend_colors.append(color_mapping.get(ct, (0.5, 0.5, 0.5)))
    
    # 创建图例元素
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=color, markersize=8, 
                                 label=ct) for ct, color in zip(unique_celltypes, legend_colors)]
    
    # 添加图例
    fig.legend(handles=legend_elements, bbox_to_anchor=(0.5, 0.02), 
              loc='upper center', ncol=4, frameon=False)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # 为图例留出空间
    
    # 保存图像
    output_file = os.path.join(output_dir, "umap_composite.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"复合UMAP图已保存到: {output_file}")


def main():
    """主函数"""
    input_dir = "/home/duxuyan/Projects/0_HiRES/output"
    output_dir = "/home/duxuyan/Projects/0_HiRES/output/graph"
    color_mapping_file = "/home/duxuyan/Projects/0_HiRES/output/color_mapping.json"
    
    # 加载所有Stage的h5ad文件
    stage_data = load_h5ad_files(input_dir)
    
    # 创建统一的颜色映射
    color_mapping = create_unified_color_mapping(stage_data, color_mapping_file)
    
    # 为所有Stage绘制UMAP图
    plot_umap_all_stages(stage_data, output_dir, color_mapping)
    
    # 生成复合UMAP图
    plot_umap_composite(stage_data, output_dir, color_mapping)
    
    print("所有Stage的UMAP可视化完成!")


if __name__ == "__main__":
    main()