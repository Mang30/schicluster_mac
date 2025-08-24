#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hi-C特化的UMAP可视化模块
基于color_mapping.json生成一致性颜色的UMAP图表
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import time

# 设置matplotlib中文字体支持
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class HiCUMAPVisualizer:
    """Hi-C特化的UMAP可视化器"""
    
    def __init__(self, color_mapping_file: str):
        """
        初始化可视化器
        
        Parameters:
        -----------
        color_mapping_file : str
            颜色映射JSON文件路径
        """
        self.color_mapping_file = color_mapping_file
        self.color_mapping = self._load_color_mapping()
        
        # 设置scanpy参数
        sc.settings.verbosity = 1  # 减少输出
        sc.settings.set_figure_params(dpi=300, facecolor='white', format='png')
        
        logger.info(f"初始化完成，加载 {len(self.color_mapping)} 种细胞类型颜色映射")
        
    def _load_color_mapping(self) -> Dict[str, Tuple[float, float, float]]:
        """加载颜色映射文件"""
        try:
            with open(self.color_mapping_file, 'r', encoding='utf-8') as f:
                color_data = json.load(f)
            
            # 转换为元组格式
            color_mapping = {}
            for celltype, rgb in color_data.items():
                color_mapping[celltype] = tuple(rgb)
                
            logger.info(f"加载颜色映射: {len(color_mapping)} 种细胞类型")
            return color_mapping
            
        except Exception as e:
            logger.error(f"加载颜色映射失败: {e}")
            raise
            
    def compute_umap_for_stage(self, adata: sc.AnnData, stage_name: str,
                              n_neighbors: int = 15, n_components: int = 2,
                              min_dist: float = 0.1, spread: float = 1.0,
                              random_state: int = 42) -> sc.AnnData:
        """
        为单个stage计算UMAP嵌入
        
        Parameters:
        -----------
        adata : sc.AnnData
            输入的AnnData对象
        stage_name : str
            stage名称
        n_neighbors : int
            UMAP邻居数
        n_components : int
            UMAP降维后维数
        min_dist : float
            UMAP最小距离参数
        spread : float
            UMAP扩散参数
        random_state : int
            随机种子
            
        Returns:
        --------
        sc.AnnData : 包含UMAP结果的AnnData对象
        """
        logger.info(f"开始计算Stage {stage_name} 的UMAP嵌入")
        start_time = time.time()
        
        try:
            # 复制数据避免修改原始对象
            adata_copy = adata.copy()
            
            # Hi-C数据特化预处理
            logger.info("进行Hi-C特化的数据预处理")
            
            # 1. 特征选择（选择高变异的bin interactions）
            if adata_copy.n_vars > 10000:  # 对于大量特征进行选择
                # 计算特征的变异系数
                from scipy.stats import variation
                feature_var = np.array([variation(adata_copy.X[:, i]) for i in range(adata_copy.n_vars)])
                
                # 选择top变异特征
                n_top_features = min(5000, adata_copy.n_vars // 2)
                top_features_idx = np.argsort(feature_var)[-n_top_features:]
                
                adata_copy = adata_copy[:, top_features_idx]
                logger.info(f"选择了 {n_top_features} 个高变异特征")
            
            # 2. PCA降维（Hi-C数据通常需要更多主成分）
            n_pcs = min(50, adata_copy.n_obs - 1, adata_copy.n_vars)
            sc.tl.pca(adata_copy, n_comps=n_pcs, random_state=random_state)
            logger.info(f"PCA降维完成: {n_pcs} 个主成分")
            
            # 3. 构建邻域图
            sc.pp.neighbors(adata_copy, 
                           n_neighbors=min(n_neighbors, adata_copy.n_obs - 1),
                           n_pcs=n_pcs, 
                           random_state=random_state)
            
            # 4. 计算UMAP
            sc.tl.umap(adata_copy, 
                      n_components=n_components,
                      min_dist=min_dist,
                      spread=spread,
                      random_state=random_state)
            
            elapsed_time = time.time() - start_time
            logger.info(f"Stage {stage_name} UMAP计算完成，耗时: {elapsed_time:.2f}秒")
            
            return adata_copy
            
        except Exception as e:
            logger.error(f"Stage {stage_name} UMAP计算失败: {e}")
            raise
            
    def plot_single_stage_umap(self, adata: sc.AnnData, stage_name: str,
                              output_path: str, figsize: Tuple[int, int] = (8, 6),
                              point_size: int = 20, alpha: float = 0.7) -> bool:
        """
        绘制单个stage的UMAP图
        
        Parameters:
        -----------
        adata : sc.AnnData
            包含UMAP结果的AnnData对象
        stage_name : str
            stage名称
        output_path : str
            输出图片路径
        figsize : Tuple[int, int]
            图片大小
        point_size : int
            散点大小
        alpha : float
            透明度
            
        Returns:
        --------
        bool : 是否成功生成图片
        """
        logger.info(f"绘制Stage {stage_name} 的UMAP图")
        
        try:
            # 创建输出目录
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            
            # 创建图形
            fig, ax = plt.subplots(figsize=figsize)
            
            # 获取细胞类型和坐标
            if 'celltype' not in adata.obs.columns:
                logger.error("AnnData对象缺少celltype信息")
                return False
                
            celltypes = adata.obs['celltype']
            umap_coords = adata.obsm['X_umap']
            
            # 按细胞类型分组绘制
            unique_celltypes = celltypes.unique()
            legend_elements = []
            
            for celltype in unique_celltypes:
                # 获取该细胞类型的索引
                mask = celltypes == celltype
                if not mask.any():
                    continue
                    
                # 获取颜色
                color = self.color_mapping.get(celltype, (0.5, 0.5, 0.5))  # 默认灰色
                
                # 绘制散点
                scatter = ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                                   c=[color], s=point_size, alpha=alpha, 
                                   label=celltype, edgecolors='none')
                
                # 添加到图例
                legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                                markerfacecolor=color, markersize=8,
                                                label=celltype))
            
            # 设置图形属性
            ax.set_xlabel('UMAP 1', fontsize=12)
            ax.set_ylabel('UMAP 2', fontsize=12)
            ax.set_title(f'Stage {stage_name}', fontsize=14, fontweight='bold')
            
            # 添加图例
            if legend_elements:
                ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), 
                         loc='upper left', frameon=False, fontsize=10)
            
            # 设置坐标轴
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.grid(True, alpha=0.3)
            
            # 保存图片
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            logger.info(f"Stage {stage_name} UMAP图已保存: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"绘制Stage {stage_name} UMAP图失败: {e}")
            return False
            
    def plot_composite_umap(self, stage_adatas: Dict[str, sc.AnnData], 
                           output_path: str, figsize: Tuple[int, int] = (20, 12),
                           point_size: int = 8, alpha: float = 0.7) -> bool:
        """
        绘制7个stage的复合UMAP图，按照FigC.png样式
        
        Parameters:
        -----------
        stage_adatas : Dict[str, sc.AnnData]
            包含各stage UMAP结果的字典
        output_path : str
            输出图片路径
        figsize : Tuple[int, int]
            图片大小
        point_size : int
            散点大小
        alpha : float
            透明度
            
        Returns:
        --------
        bool : 是否成功生成复合图
        """
        logger.info("绘制7-stage复合UMAP图")
        
        try:
            # 定义stage顺序（按照FigC.png）
            stage_order = ['E70', 'E75', 'E80', 'E85', 'E95', 'EX05', 'EX15']
            # 实际存在的stage
            available_stages = [s for s in stage_order if s in stage_adatas]
            
            if not available_stages:
                logger.error("没有可用的stage数据")
                return False
                
            logger.info(f"绘制 {len(available_stages)} 个stage: {available_stages}")
            
            # 创建子图布局 (2行4列，最后一个空位)
            fig, axes = plt.subplots(2, 4, figsize=figsize)
            axes = axes.flatten()
            
            # 隐藏最后一个子图
            if len(available_stages) < 8:
                axes[-1].set_visible(False)
                
            # 收集所有细胞类型用于统一图例
            all_celltypes = set()
            for adata in stage_adatas.values():
                if 'celltype' in adata.obs.columns:
                    all_celltypes.update(adata.obs['celltype'].unique())
            all_celltypes = sorted(list(all_celltypes))
            
            # 为每个stage绘制子图
            for i, stage_name in enumerate(available_stages):
                ax = axes[i]
                adata = stage_adatas[stage_name]
                
                if 'X_umap' not in adata.obsm:
                    logger.warning(f"Stage {stage_name} 缺少UMAP结果")
                    ax.text(0.5, 0.5, f'{stage_name}\n(No UMAP)', 
                           ha='center', va='center', transform=ax.transAxes)
                    ax.set_title(stage_name, fontsize=12, fontweight='bold')
                    continue
                
                # 获取数据
                celltypes = adata.obs.get('celltype', pd.Series(['Unknown'] * adata.n_obs))
                umap_coords = adata.obsm['X_umap']
                
                # 按细胞类型绘制
                for celltype in celltypes.unique():
                    mask = celltypes == celltype
                    if not mask.any():
                        continue
                        
                    color = self.color_mapping.get(celltype, (0.5, 0.5, 0.5))
                    
                    ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                             c=[color], s=point_size, alpha=alpha, 
                             edgecolors='none')
                
                # 设置子图属性
                ax.set_title(stage_name, fontsize=12, fontweight='bold')
                ax.set_xlabel('Hi-C UMAP 1', fontsize=10)
                ax.set_ylabel('Hi-C UMAP 2', fontsize=10) if i % 4 == 0 else ax.set_ylabel('')
                
                # 移除边框
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.tick_params(labelsize=8)
                
            # 创建统一图例
            legend_elements = []
            for celltype in all_celltypes:
                if celltype in self.color_mapping:
                    color = self.color_mapping[celltype]
                    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                                    markerfacecolor=color, markersize=8,
                                                    label=celltype))
            
            # 添加图例到图片底部
            if legend_elements:
                fig.legend(handles=legend_elements, 
                          loc='lower center', 
                          bbox_to_anchor=(0.5, -0.02),
                          ncol=min(6, len(legend_elements)), 
                          frameon=False, 
                          fontsize=10)
            
            # 调整布局
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.15)  # 为图例留出空间
            
            # 保存图片
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            logger.info(f"复合UMAP图已保存: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"绘制复合UMAP图失败: {e}")
            return False
            
    def process_all_stages(self, h5ad_dir: str, output_dir: str,
                          stages: Optional[List[str]] = None) -> Dict[str, bool]:
        """
        处理所有stage的UMAP可视化
        
        Parameters:
        -----------
        h5ad_dir : str
            h5ad文件目录
        output_dir : str
            输出目录
        stages : List[str], optional
            要处理的stage列表
            
        Returns:
        --------
        Dict[str, bool] : 各stage处理结果
        """
        logger.info("开始处理所有stage的UMAP可视化")
        
        # 默认stage列表
        if stages is None:
            stages = ['E70', 'E75', 'E80', 'E85', 'E95', 'EX05', 'EX15']
            
        # 创建输出目录
        single_plots_dir = os.path.join(output_dir, 'single_stage_umaps')
        os.makedirs(single_plots_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        stage_results = {}
        stage_adatas = {}
        
        # 处理每个stage
        for stage in stages:
            h5ad_file = os.path.join(h5ad_dir, f"{stage}_hic.h5ad")
            
            if not os.path.exists(h5ad_file):
                logger.warning(f"h5ad文件不存在: {h5ad_file}")
                stage_results[stage] = False
                continue
                
            try:
                # 加载数据
                logger.info(f"处理Stage {stage}")
                adata = sc.read_h5ad(h5ad_file)
                
                # 检查细胞数量
                if adata.n_obs < 3:
                    logger.warning(f"Stage {stage} 细胞数量不足: {adata.n_obs}")
                    stage_results[stage] = False
                    continue
                
                # 计算UMAP
                adata_with_umap = self.compute_umap_for_stage(adata, stage)
                
                # 绘制单独的图
                single_output_path = os.path.join(single_plots_dir, f"umap_{stage}.png")
                single_success = self.plot_single_stage_umap(adata_with_umap, stage, single_output_path)
                
                if single_success:
                    stage_adatas[stage] = adata_with_umap
                    stage_results[stage] = True
                    logger.info(f"Stage {stage} 处理成功")
                else:
                    stage_results[stage] = False
                    logger.error(f"Stage {stage} 单图绘制失败")
                    
            except Exception as e:
                logger.error(f"处理Stage {stage} 失败: {e}")
                stage_results[stage] = False
                
        # 绘制复合图
        if stage_adatas:
            composite_output_path = os.path.join(output_dir, "composite_umap_7stages.png")
            composite_success = self.plot_composite_umap(stage_adatas, composite_output_path)
            
            if composite_success:
                logger.info("复合UMAP图生成成功")
            else:
                logger.error("复合UMAP图生成失败")
                
        logger.info(f"所有stage处理完成: {sum(stage_results.values())}/{len(stages)} 成功")
        return stage_results


def main():
    """测试主函数"""
    # 配置路径
    color_mapping_file = "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/plotumap/color_mapping.json"
    h5ad_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/hic_h5ad_files"
    output_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/umap_plots"
    
    # 创建可视化器
    visualizer = HiCUMAPVisualizer(color_mapping_file)
    
    # 测试阶段
    test_stages = ['E70', 'E80', 'EX05']
    
    # 处理所有stage
    results = visualizer.process_all_stages(h5ad_dir, output_dir, test_stages)
    
    logger.info("UMAP可视化测试完成")
    for stage, success in results.items():
        logger.info(f"Stage {stage}: {'成功' if success else '失败'}")


if __name__ == "__main__":
    main()