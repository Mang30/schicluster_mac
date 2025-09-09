#!/usr/bin/env python3
"""
plot_stage_E75_umap.py - 专门用于绘制stage_E75.h5ad的UMAP图
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def load_color_mapping(color_mapping_file):
    """
    加载颜色映射文件
    """
    try:
        with open(color_mapping_file, 'r') as f:
            color_mapping = json.load(f)
        print(f"成功加载颜色映射，包含 {len(color_mapping)} 个细胞类型")
        return color_mapping
    except Exception as e:
        print(f"加载颜色映射文件失败: {e}")
        return {}

def plot_stage_E75_umap():
    """
    专门为stage_E75.h5ad文件绘制UMAP图
    """
    # 设置文件路径
    h5ad_file = "/Volumes/SumSung500/CSU/0_HiRES/4_contact_decay_profile_2_h5ad/json2h5ad/output/stage_E75.h5ad"
    color_mapping_file = "/Volumes/SumSung500/CSU/0_HiRES/5_plot_umap/color_mapping.json"
    
    # 创建stage-specific输出目录
    base_output_dir = "/Volumes/SumSung500/CSU/0_HiRES/5_plot_umap/output"
    stage_name = "stage_E75"
    output_dir = os.path.join(base_output_dir, stage_name)
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    print(f"输出目录: {output_dir}")
    
    # 检查H5AD文件是否存在
    if not os.path.exists(h5ad_file):
        print(f"错误: H5AD文件不存在: {h5ad_file}")
        return
    
    # 加载颜色映射
    color_mapping = load_color_mapping(color_mapping_file)
    
    print(f"正在读取H5AD文件: {h5ad_file}")
    
    try:
        # 读取h5ad文件
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据形状: {adata.shape}")
        print(f"观测值列: {list(adata.obs.columns)}")
        print(f"变量列: {list(adata.var.columns) if adata.var.columns.size > 0 else '无'}")
        
        # 显示细胞类型分布
        if 'celltype' in adata.obs.columns:
            print(f"\\n细胞类型分布:")
            print(adata.obs['celltype'].value_counts())
        elif 'Celltype' in adata.obs.columns:
            print(f"\\n细胞类型分布:")
            print(adata.obs['Celltype'].value_counts())
        
        # 设置scanpy参数
        sc.settings.verbosity = 1
        sc.settings.set_figure_params(dpi=80, facecolor='white')
        
        # 检查是否需要数据预处理
        print("检查数据质量...")
        print(f"数据矩阵形状: {adata.X.shape}")
        print(f"数据类型: {adata.X.dtype}")
        print(f"数据范围: {adata.X.min():.6f} - {adata.X.max():.6f}")
        
        # 检查NaN值
        nan_count = np.isnan(adata.X).sum()
        print(f"NaN值数量: {nan_count}")
        
        if nan_count > 0:
            print(f"发现 {nan_count} 个NaN值，正在处理...")
            # 用0替换NaN值（对于contact decay profile，0是合理的默认值）
            adata.X = np.nan_to_num(adata.X, nan=0.0)
            print("NaN值已替换为0")
        
        # 检查是否有无限值
        inf_count = np.isinf(adata.X).sum()
        if inf_count > 0:
            print(f"发现 {inf_count} 个无限值，正在处理...")
            adata.X = np.nan_to_num(adata.X, posinf=adata.X[~np.isinf(adata.X)].max(), 
                                              neginf=adata.X[~np.isinf(adata.X)].min())
            print("无限值已处理")
        
        if adata.X.max() > 100:  # 如果数据看起来是原始counts
            print("正在进行数据标准化...")
            # 保存原始数据
            adata.raw = adata
            # 标准化到每个细胞总counts为1e4
            sc.pp.normalize_total(adata, target_sum=1e4)
            # log转换
            sc.pp.log1p(adata)
        else:
            print("数据已经标准化，跳过标准化步骤")
        
        # 再次检查处理后的数据
        print(f"处理后数据范围: {adata.X.min():.6f} - {adata.X.max():.6f}")
        print(f"处理后NaN数量: {np.isnan(adata.X).sum()}")
        print(f"处理后Inf数量: {np.isinf(adata.X).sum()}")
        
        # 计算PCA
        print("计算主成分分析 (PCA)...")
        n_features = adata.shape[1]
        n_samples = adata.shape[0]
        n_comps = min(50, min(n_samples, n_features) - 1)
        
        if n_comps < 2:
            print("警告: 特征数量太少，无法进行PCA分析")
            return
            
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
        
        # 计算邻居图
        print("计算邻居图...")
        n_neighbors = min(15, n_samples - 1)
        if n_neighbors < 2:
            n_neighbors = 2
        
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(40, n_comps))
        
        # 计算UMAP
        print("计算UMAP嵌入...")
        sc.tl.umap(adata)
        
        # 创建图形
        plt.style.use('default')
        
        # 1. 按细胞类型绘制UMAP
        celltype_col = None
        if 'celltype' in adata.obs.columns:
            celltype_col = 'celltype'
        elif 'Celltype' in adata.obs.columns:
            celltype_col = 'Celltype'
        
        if celltype_col:
            print(f"绘制按{celltype_col}着色的UMAP图...")
            
            # 获取唯一的细胞类型
            unique_celltypes = adata.obs[celltype_col].dropna().unique()
            print(f"发现 {len(unique_celltypes)} 种细胞类型: {list(unique_celltypes)}")
            
            # 创建图形
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # 获取UMAP坐标
            umap_coords = adata.obsm['X_umap']
            
            # 为每种细胞类型绘制散点
            for i, celltype in enumerate(unique_celltypes):
                mask = adata.obs[celltype_col] == celltype
                
                # 获取颜色
                if celltype in color_mapping:
                    color = color_mapping[celltype]
                else:
                    # 使用默认颜色
                    color = plt.cm.tab20(i % 20)
                
                # 绘制散点
                ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1], 
                          c=[color], label=celltype, s=60, alpha=0.8, edgecolors='none')
            
            ax.set_xlabel('UMAP1', fontsize=12)
            ax.set_ylabel('UMAP2', fontsize=12)
            ax.set_title(f'Stage E75 - UMAP by Cell Type\\n({adata.shape[0]} cells)', fontsize=14)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # 保存图像
            output_file = os.path.join(output_dir, "umap_celltype.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            print(f"保存UMAP图 (按细胞类型): {output_file}")
        
        # 2. 如果有发育阶段信息，绘制按stage着色的图
        stage_col = None
        if 'stage' in adata.obs.columns:
            stage_col = 'stage'
        elif 'Stage' in adata.obs.columns:
            stage_col = 'Stage'
        
        if stage_col:
            print(f"绘制按{stage_col}着色的UMAP图...")
            
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # 使用scanpy的绘图函数
            sc.pl.umap(adata, color=stage_col, ax=ax, show=False, legend_loc='right margin')
            plt.title(f'Stage E75 - UMAP by Developmental Stage\\n({adata.shape[0]} cells)', fontsize=14)
            plt.tight_layout()
            
            # 保存图像
            output_file = os.path.join(output_dir, "umap_stage.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            print(f"保存UMAP图 (按发育阶段): {output_file}")
        
        # 3. 如果有细胞周期信息，绘制按细胞周期着色的图
        cellcycle_col = None
        if 'cellcycle_phase' in adata.obs.columns:
            cellcycle_col = 'cellcycle_phase'
        elif 'Cellcycle phase' in adata.obs.columns:
            cellcycle_col = 'Cellcycle phase'
        
        if cellcycle_col:
            print(f"绘制按{cellcycle_col}着色的UMAP图...")
            
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # 使用scanpy的绘图函数
            sc.pl.umap(adata, color=cellcycle_col, ax=ax, show=False, legend_loc='right margin')
            plt.title(f'Stage E75 - UMAP by Cell Cycle Phase\\n({adata.shape[0]} cells)', fontsize=14)
            plt.tight_layout()
            
            # 保存图像
            output_file = os.path.join(output_dir, "umap_cellcycle.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            print(f"保存UMAP图 (按细胞周期): {output_file}")
        
        # 4. 绘制一些数值型特征
        numeric_features = []
        for col in adata.obs.columns:
            if col in ['total_contacts', 'power_law_slope', 'g1s_score', 'g2m_score', 'repli_score']:
                if col in adata.obs.columns and adata.obs[col].dtype in ['float64', 'int64']:
                    numeric_features.append(col)
        
        if numeric_features:
            print(f"绘制数值特征的UMAP图: {numeric_features}")
            
            # 创建多个子图
            n_features = len(numeric_features)
            n_cols = min(3, n_features)
            n_rows = (n_features + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
            if n_features == 1:
                axes = [axes]
            elif n_rows == 1:
                axes = axes
            else:
                axes = axes.flatten()
            
            for i, feature in enumerate(numeric_features):
                ax = axes[i] if n_features > 1 else axes[0]
                sc.pl.umap(adata, color=feature, ax=ax, show=False, colorbar_loc='right')
                ax.set_title(f'{feature}')
            
            # 隐藏多余的子图
            for i in range(n_features, len(axes)):
                axes[i].set_visible(False)
            
            plt.suptitle(f'Stage E75 - UMAP by Numeric Features', fontsize=16)
            plt.tight_layout()
            
            # 保存图像
            output_file = os.path.join(output_dir, "umap_numeric_features.png")
            plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            print(f"保存UMAP图 (数值特征): {output_file}")
        
        print(f"\\n所有UMAP图已保存到目录: {output_dir}")
        print("绘图完成!")
        
        # 保存数据摘要
        summary_file = os.path.join(output_dir, "data_summary.txt")
        with open(summary_file, 'w') as f:
            f.write(f"Stage E75 数据摘要\\n")
            f.write(f"=" * 50 + "\\n")
            f.write(f"细胞数量: {adata.shape[0]}\\n")
            f.write(f"特征数量: {adata.shape[1]}\\n")
            f.write(f"观测值列: {list(adata.obs.columns)}\\n")
            f.write(f"\\n细胞类型分布:\\n")
            if celltype_col:
                for ct, count in adata.obs[celltype_col].value_counts().items():
                    f.write(f"  {ct}: {count}\\n")
        
        print(f"数据摘要已保存: {summary_file}")
        
    except Exception as e:
        print(f"处理H5AD文件时出错: {e}")
        import traceback
        traceback.print_exc()

def main():
    """
    主函数
    """
    print("开始绘制Stage E75的UMAP图...")
    plot_stage_E75_umap()

if __name__ == "__main__":
    main()
