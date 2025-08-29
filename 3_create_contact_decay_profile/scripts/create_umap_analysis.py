#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
接触衰减特征UMAP分析脚本

功能：
1. 将长格式的衰减曲线数据转换为细胞×特征矩阵
2. 进行特征预处理和降维
3. 生成UMAP可视化图表
4. 保存可用于下游分析的特征矩阵

作者：Claude Code Assistant  
日期：2025-08-29
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
import warnings
warnings.filterwarnings('ignore')

try:
    import umap
    from sklearn.preprocessing import StandardScaler, RobustScaler
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
except ImportError as e:
    print(f"错误：缺少必要的包: {e}")
    print("请安装：pip install umap-learn scikit-learn")
    sys.exit(1)

# 配置matplotlib中文显示
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class ContactDecayUMAPAnalyzer:
    """接触衰减特征UMAP分析器"""
    
    def __init__(self, data_file: str, output_dir: str):
        """
        初始化UMAP分析器
        
        Args:
            data_file: 长格式的衰减曲线数据文件路径
            output_dir: 输出目录
        """
        self.data_file = Path(data_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置日志
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # 数据存储
        self.raw_data = None
        self.feature_matrix = None
        self.cell_metadata = None
        self.umap_embedding = None
        self.pca_embedding = None
        
    def load_decay_data(self) -> bool:
        """加载衰减曲线数据"""
        try:
            if not self.data_file.exists():
                raise FileNotFoundError(f"数据文件不存在: {self.data_file}")
            
            self.logger.info(f"加载数据: {self.data_file}")
            self.raw_data = pd.read_csv(self.data_file)
            
            self.logger.info(f"数据形状: {self.raw_data.shape}")
            self.logger.info(f"发现细胞数: {self.raw_data.groupby(['stage', 'sample']).ngroups}")
            self.logger.info(f"发现发育阶段: {sorted(self.raw_data['stage'].unique())}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"加载数据失败: {e}")
            return False
    
    def create_feature_matrix(self, distance_range: Optional[Tuple[float, float]] = None,
                             log_transform: bool = True, 
                             feature_selection: str = 'all') -> bool:
        """
        创建细胞×特征矩阵
        
        Args:
            distance_range: 距离范围过滤 (min_kb, max_kb)
            log_transform: 是否使用对数变换的接触频率
            feature_selection: 特征选择策略 ('all', 'uniform', 'log_uniform')
            
        Returns:
            bool: 是否成功创建
        """
        try:
            if self.raw_data is None:
                raise ValueError("请先加载数据")
            
            self.logger.info("创建细胞×特征矩阵...")
            
            # 数据预处理
            data = self.raw_data.copy()
            
            # 距离范围过滤
            if distance_range:
                min_kb, max_kb = distance_range
                data = data[
                    (data['distance_kb'] >= min_kb) & 
                    (data['distance_kb'] <= max_kb)
                ]
                self.logger.info(f"距离过滤: {min_kb}-{max_kb} kb, 剩余数据: {len(data)}")
            
            # 选择接触频率特征
            if log_transform:
                # 使用对数变换的接触频率，处理无穷大值
                data['feature_value'] = np.log10(data['contact_frequency'] + 1e-10)
                data['feature_value'] = np.nan_to_num(data['feature_value'], 
                                                     nan=0.0, posinf=0.0, neginf=-10.0)
                self.logger.info("使用对数变换的接触频率")
            else:
                data['feature_value'] = data['contact_frequency']
                self.logger.info("使用原始接触频率")
            
            # 特征选择
            if feature_selection == 'uniform':
                # 均匀采样距离点
                unique_distances = sorted(data['distance_kb'].unique())
                n_features = min(200, len(unique_distances))  # 最多200个特征
                selected_distances = unique_distances[::len(unique_distances)//n_features]
                data = data[data['distance_kb'].isin(selected_distances)]
                self.logger.info(f"均匀采样: 选择了 {len(selected_distances)} 个距离点")
                
            elif feature_selection == 'log_uniform':
                # 对数均匀采样
                unique_distances = sorted(data['distance_kb'].unique())
                log_distances = np.log10(unique_distances)
                n_features = min(200, len(unique_distances))
                selected_log = np.linspace(log_distances[0], log_distances[-1], n_features)
                
                # 找到最接近的实际距离
                selected_distances = []
                for target_log in selected_log:
                    target_dist = 10**target_log
                    closest_dist = min(unique_distances, key=lambda x: abs(x - target_dist))
                    if closest_dist not in selected_distances:
                        selected_distances.append(closest_dist)
                
                data = data[data['distance_kb'].isin(selected_distances)]
                self.logger.info(f"对数均匀采样: 选择了 {len(selected_distances)} 个距离点")
            
            # 创建透视表（细胞×特征矩阵）
            pivot_table = data.pivot_table(
                index=['stage', 'sample'],
                columns='distance_kb', 
                values='feature_value',
                aggfunc='first'  # 每个细胞每个距离只有一个值
            )
            
            # 处理缺失值
            pivot_table = pivot_table.fillna(0)
            
            # 提取特征矩阵和元数据
            self.feature_matrix = pivot_table.values
            self.cell_metadata = pd.DataFrame({
                'stage': [idx[0] for idx in pivot_table.index],
                'sample': [idx[1] for idx in pivot_table.index],
                'cell_id': [f"{idx[0]}_{idx[1]}" for idx in pivot_table.index]
            })
            
            self.logger.info(f"特征矩阵形状: {self.feature_matrix.shape}")
            self.logger.info(f"细胞数量: {len(self.cell_metadata)}")
            self.logger.info(f"特征维度: {self.feature_matrix.shape[1]}")
            
            # 保存特征矩阵
            feature_df = pd.DataFrame(
                self.feature_matrix, 
                columns=[f"distance_{int(col)}kb" for col in pivot_table.columns],
                index=self.cell_metadata['cell_id']
            )
            feature_df = pd.concat([self.cell_metadata.set_index('cell_id'), feature_df], axis=1)
            
            output_file = self.output_dir / 'cell_feature_matrix.csv'
            feature_df.to_csv(output_file)
            self.logger.info(f"特征矩阵已保存: {output_file}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"创建特征矩阵失败: {e}")
            return False
    
    def perform_dimensionality_reduction(self, method: str = 'umap', 
                                       n_components: int = 2,
                                       preprocessing: str = 'standard') -> bool:
        """
        进行降维分析
        
        Args:
            method: 降维方法 ('umap', 'pca', 'tsne')
            n_components: 降维后的维度
            preprocessing: 预处理方法 ('standard', 'robust', 'none')
            
        Returns:
            bool: 是否成功
        """
        try:
            if self.feature_matrix is None:
                raise ValueError("请先创建特征矩阵")
            
            self.logger.info(f"开始降维分析: {method}")
            
            # 数据预处理
            if preprocessing == 'standard':
                scaler = StandardScaler()
                X_scaled = scaler.fit_transform(self.feature_matrix)
                self.logger.info("使用标准化预处理")
            elif preprocessing == 'robust':
                scaler = RobustScaler()
                X_scaled = scaler.fit_transform(self.feature_matrix)
                self.logger.info("使用鲁棒标准化预处理")
            else:
                X_scaled = self.feature_matrix
                self.logger.info("不进行预处理")
            
            # 降维
            if method.lower() == 'umap':
                reducer = umap.UMAP(
                    n_components=n_components,
                    n_neighbors=15,
                    min_dist=0.1,
                    metric='euclidean',
                    random_state=42
                )
                embedding = reducer.fit_transform(X_scaled)
                self.umap_embedding = embedding
                self.logger.info(f"UMAP降维完成: {self.feature_matrix.shape} -> {embedding.shape}")
                
            elif method.lower() == 'pca':
                reducer = PCA(n_components=n_components, random_state=42)
                embedding = reducer.fit_transform(X_scaled)
                self.pca_embedding = embedding
                
                # 保存解释方差比例
                explained_var = reducer.explained_variance_ratio_
                self.logger.info(f"PCA降维完成: 前{n_components}个主成分解释了 {explained_var.sum():.3f} 的方差")
                
            elif method.lower() == 'tsne':
                if X_scaled.shape[1] > 50:  # t-SNE对高维数据先进行PCA预处理
                    pca = PCA(n_components=50, random_state=42)
                    X_scaled = pca.fit_transform(X_scaled)
                    self.logger.info("t-SNE前进行PCA预处理至50维")
                
                reducer = TSNE(n_components=n_components, random_state=42, 
                             perplexity=min(30, len(X_scaled)//4))
                embedding = reducer.fit_transform(X_scaled)
                self.logger.info(f"t-SNE降维完成")
            
            # 保存降维结果
            embedding_df = pd.DataFrame(
                embedding,
                columns=[f'{method.upper()}{i+1}' for i in range(n_components)],
                index=self.cell_metadata['cell_id']
            )
            embedding_df = pd.concat([self.cell_metadata.set_index('cell_id'), embedding_df], axis=1)
            
            output_file = self.output_dir / f'{method}_embedding.csv'
            embedding_df.to_csv(output_file)
            self.logger.info(f"降维结果已保存: {output_file}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"降维分析失败: {e}")
            return False
    
    def create_umap_plots(self, color_by: str = 'stage', save_plots: bool = True) -> List[str]:
        """
        创建UMAP可视化图表
        
        Args:
            color_by: 着色依据 ('stage', 'sample')
            save_plots: 是否保存图表
            
        Returns:
            List[str]: 生成的图表文件路径
        """
        try:
            if self.umap_embedding is None:
                raise ValueError("请先进行UMAP降维")
            
            plot_files = []
            
            # 准备数据
            plot_data = pd.DataFrame({
                'UMAP1': self.umap_embedding[:, 0],
                'UMAP2': self.umap_embedding[:, 1],
                'stage': self.cell_metadata['stage'],
                'sample': self.cell_metadata['sample'],
                'cell_id': self.cell_metadata['cell_id']
            })
            
            # 设置颜色方案
            if color_by == 'stage':
                unique_stages = sorted(plot_data['stage'].unique())
                colors = plt.cm.Set3(np.linspace(0, 1, len(unique_stages)))
                color_map = dict(zip(unique_stages, colors))
                
            # 1. 基本UMAP图
            fig, ax = plt.subplots(figsize=(10, 8))
            
            for stage in sorted(plot_data['stage'].unique()):
                stage_data = plot_data[plot_data['stage'] == stage]
                ax.scatter(stage_data['UMAP1'], stage_data['UMAP2'],
                          c=[color_map[stage]], label=f'阶段 {stage}',
                          alpha=0.7, s=50)
            
            ax.set_xlabel('UMAP1', fontsize=12)
            ax.set_ylabel('UMAP2', fontsize=12)
            ax.set_title('Hi-C接触衰减特征 UMAP 分析', fontsize=14)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            if save_plots:
                plot_file = self.output_dir / 'umap_by_stage.png'
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                plot_files.append(str(plot_file))
                
            plt.show()
            
            # 2. 密度图
            fig, axes = plt.subplots(2, 2, figsize=(16, 12))
            
            # 整体密度
            axes[0,0].hexbin(plot_data['UMAP1'], plot_data['UMAP2'], gridsize=30, cmap='Blues')
            axes[0,0].set_title('细胞分布密度')
            axes[0,0].set_xlabel('UMAP1')
            axes[0,0].set_ylabel('UMAP2')
            
            # 各阶段分布
            for i, stage in enumerate(sorted(plot_data['stage'].unique())[:3]):
                ax = axes[0,1] if i == 0 else axes[1, i-1]
                stage_data = plot_data[plot_data['stage'] == stage]
                ax.scatter(stage_data['UMAP1'], stage_data['UMAP2'],
                          c=color_map[stage], alpha=0.6, s=30)
                ax.set_title(f'阶段 {stage}')
                ax.set_xlabel('UMAP1')
                ax.set_ylabel('UMAP2')
                ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            if save_plots:
                plot_file = self.output_dir / 'umap_density_analysis.png'
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                plot_files.append(str(plot_file))
                
            plt.show()
            
            # 3. 统计摘要
            stage_stats = plot_data.groupby('stage').agg({
                'UMAP1': ['mean', 'std'],
                'UMAP2': ['mean', 'std']
            }).round(3)
            
            print("\n📊 各发育阶段UMAP坐标统计:")
            print(stage_stats)
            
            # 保存统计结果
            stats_file = self.output_dir / 'umap_statistics.csv' 
            stage_stats.to_csv(stats_file)
            
            self.logger.info(f"生成了 {len(plot_files)} 个UMAP图表")
            return plot_files
            
        except Exception as e:
            self.logger.error(f"创建UMAP图表失败: {e}")
            return []
    
    def run_complete_analysis(self, distance_range: Tuple[float, float] = (100, 10000),
                            feature_selection: str = 'log_uniform') -> Dict:
        """
        运行完整的UMAP分析流程
        
        Args:
            distance_range: 分析的距离范围 (kb)
            feature_selection: 特征选择策略
            
        Returns:
            Dict: 分析结果摘要
        """
        try:
            self.logger.info("开始完整UMAP分析流程...")
            
            # 1. 加载数据
            if not self.load_decay_data():
                return {'success': False, 'error': '数据加载失败'}
            
            # 2. 创建特征矩阵
            if not self.create_feature_matrix(
                distance_range=distance_range,
                log_transform=True,
                feature_selection=feature_selection
            ):
                return {'success': False, 'error': '特征矩阵创建失败'}
            
            # 3. 执行UMAP降维
            if not self.perform_dimensionality_reduction(method='umap'):
                return {'success': False, 'error': 'UMAP降维失败'}
            
            # 4. 生成可视化图表
            plot_files = self.create_umap_plots()
            
            # 5. 返回结果摘要
            result_summary = {
                'success': True,
                'n_cells': len(self.cell_metadata),
                'n_features': self.feature_matrix.shape[1],
                'distance_range': distance_range,
                'feature_selection': feature_selection,
                'stages': sorted(self.cell_metadata['stage'].unique()),
                'output_dir': str(self.output_dir),
                'plot_files': plot_files
            }
            
            self.logger.info("UMAP分析流程完成!")
            return result_summary
            
        except Exception as e:
            self.logger.error(f"完整分析失败: {e}")
            return {'success': False, 'error': str(e)}


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='接触衰减特征UMAP分析工具')
    parser.add_argument('-i', '--input', required=True, 
                       help='输入的衰减曲线数据文件 (all_stages_decay_profiles.csv)')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    parser.add_argument('--min-distance', type=float, default=100,
                       help='最小分析距离 (kb)')
    parser.add_argument('--max-distance', type=float, default=10000,
                       help='最大分析距离 (kb)')
    parser.add_argument('--feature-selection', choices=['all', 'uniform', 'log_uniform'],
                       default='log_uniform', help='特征选择策略')
    
    args = parser.parse_args()
    
    # 创建分析器
    analyzer = ContactDecayUMAPAnalyzer(
        data_file=args.input,
        output_dir=args.output
    )
    
    # 运行分析
    results = analyzer.run_complete_analysis(
        distance_range=(args.min_distance, args.max_distance),
        feature_selection=args.feature_selection
    )
    
    # 输出结果
    if results['success']:
        print(f"\n🎉 UMAP分析完成！")
        print(f"📊 细胞数量: {results['n_cells']}")
        print(f"🔍 特征维度: {results['n_features']}")  
        print(f"🧬 发育阶段: {results['stages']}")
        print(f"📁 输出目录: {results['output_dir']}")
        print(f"📈 生成图表: {len(results['plot_files'])} 个")
    else:
        print(f"❌ 分析失败: {results['error']}")
        sys.exit(1)


if __name__ == "__main__":
    main()