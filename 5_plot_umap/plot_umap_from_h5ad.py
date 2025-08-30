#!/usr/bin/env python3
"""
2_plot_umap_from_h5ad.py - 从带元数据的h5ad文件绘制UMAP图
使用统一的颜色映射确保相同细胞类型在不同stage中颜色一致
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
from pathlib import Path

def load_color_mapping(color_mapping_file):
    """
    加载颜色映射
    """
    with open(color_mapping_file, 'r') as f:
        color_mapping = json.load(f)
    
    # 转换颜色值为matplotlib格式 (0-1范围)
    converted_mapping = {}
    for celltype, color in color_mapping.items():
        if isinstance(color, list) and len(color) == 3:
            # 将0-255范围转换为0-1范围
            converted_mapping[celltype] = [c for c in color]
        else:
            converted_mapping[celltype] = color
    
    return converted_mapping

def plot_umap_by_stage(input_dir, color_mapping_file, output_dir, stage=None):
    """
    为每个stage绘制UMAP图
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 加载颜色映射
    color_mapping = load_color_mapping(color_mapping_file)
    print(f"加载了 {len(color_mapping)} 个颜色映射")
    
    # 根据是否指定stage来确定要处理的文件
    if stage:
        # 处理特定stage的文件
        h5ad_files = list(Path(input_dir).glob(f"{stage}_decay_profiles.h5ad"))
        print(f"找到 {len(h5ad_files)} 个 {stage} 的h5ad文件")
    else:
        # 获取所有h5ad文件
        h5ad_files = list(Path(input_dir).glob("*.h5ad"))
        print(f"找到 {len(h5ad_files)} 个h5ad文件")
    
    # 为每个文件绘制UMAP图
    for h5ad_file in h5ad_files:
        print(f"\n处理文件: {h5ad_file}")
        
        # 读取h5ad文件
        adata = sc.read(str(h5ad_file))
        print(f"数据形状: {adata.shape}")
        
        # 检查是否有UMAP坐标，如果没有则计算
        if 'X_umap' not in adata.obsm:
            print("计算UMAP坐标...")
            # 确保数据被正确设置
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            # 计算PCA和UMAP
            # 调整n_comps和n_pcs以适应小数据集
            n_features = adata.shape[1]
            n_samples = adata.shape[0]
            n_comps = min(50, min(n_samples, n_features) - 1)
            n_pcs = min(20, n_comps)
            
            # 确保n_neighbors不超过样本数-1
            n_neighbors = min(10, n_samples - 1)
            if n_neighbors < 2:
                n_neighbors = 2
            
            print(f"使用参数: n_comps={n_comps}, n_pcs={n_pcs}, n_neighbors={n_neighbors}")
            
            sc.pp.pca(adata, n_comps=n_comps)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            sc.tl.umap(adata)
        
        # 绘制按Celltype着色的UMAP图
        if 'Celltype' in adata.obs.columns:
            # 获取存在的细胞类型
            existing_celltypes = adata.obs['Celltype'].dropna().unique()
            
            # 创建颜色列表
            colors = []
            for ct in existing_celltypes:
                if ct in color_mapping:
                    color_val = color_mapping[ct]
                    # 如果颜色值是0-1范围的浮点数列表，直接使用
                    if isinstance(color_val, list) and len(color_val) == 3:
                        colors.append(color_val)
                    else:
                        colors.append(color_val)
                else:
                    # 如果没有预定义颜色，使用默认颜色
                    colors.append(plt.cm.tab20(len(colors) % 20))
            
            # 绘制UMAP图
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # 使用自定义颜色绘制
            # 确保颜色映射正确格式化
            color_palette = dict(zip(existing_celltypes, colors))
            # 使用plot函数直接绘制，避免scanpy的颜色处理问题
            sc.tl.umap(adata)  # 确保UMAP坐标已计算
            x_coords = adata.obsm['X_umap'][:, 0]
            y_coords = adata.obsm['X_umap'][:, 1]
            celltypes = adata.obs['Celltype']
            
            # 为每个细胞类型绘制散点图
            for i, ct in enumerate(existing_celltypes):
                mask = celltypes == ct
                color = color_palette[ct] if ct in color_palette else plt.cm.tab20(i % 20)
                ax.scatter(x_coords[mask], y_coords[mask], label=ct, color=color, s=50)
            
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title(f'UMAP Plot - {h5ad_file.stem}')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            # 保存图像
            output_file = os.path.join(output_dir, f"{h5ad_file.stem}_umap_celltype.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"保存UMAP图: {output_file}")
        
        # 绘制按Stage着色的UMAP图
        if 'Stage' in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(10, 8))
            sc.pl.umap(adata, color='Stage', ax=ax, show=False)
            plt.title(f'UMAP Plot - {h5ad_file.stem} (by Stage)')
            plt.tight_layout()
            
            # 保存图像
            output_file = os.path.join(output_dir, f"{h5ad_file.stem}_umap_stage.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"保存Stage UMAP图: {output_file}")
    
    print(f"\n所有UMAP图已保存到: {output_dir}")

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='从h5ad文件绘制UMAP图')
    parser.add_argument('--input_dir', type=str, required=True, 
                        help='输入h5ad文件的目录路径')
    parser.add_argument('--color_mapping_file', type=str, required=True, 
                        help='颜色映射文件路径')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='输出UMAP图的目录路径')
    parser.add_argument('--stage', type=str, required=False,
                        help='特定stage名称 (例如: stage_E75)。如果未指定，则处理所有stage文件')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 设置路径
    input_dir = args.input_dir
    color_mapping_file = args.color_mapping_file
    output_dir = args.output_dir
    
    # 检查输入目录
    if not os.path.exists(input_dir):
        print(f"错误: 输入目录不存在: {input_dir}")
        return
    
    # 检查颜色映射文件
    if not os.path.exists(color_mapping_file):
        print(f"错误: 颜色映射文件不存在: {color_mapping_file}")
        return
    
    # 如果指定了stage，则创建新的输入目录路径
    if args.stage:
        # 对于UMAP绘图，我们假设输入目录已经包含了带元数据的文件
        # 所以我们仍然使用相同的输入目录，但在plot_umap_by_stage函数中过滤文件
        plot_umap_by_stage(input_dir, color_mapping_file, output_dir, args.stage)
    else:
        plot_umap_by_stage(input_dir, color_mapping_file, output_dir)
    
    print("UMAP绘图任务完成!")

if __name__ == "__main__":
    main()