#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
接触衰减曲线分析核心模块
用于分析Hi-C .cool格式数据的接触衰减特征

功能：
1. 从.cool文件读取Hi-C接触矩阵
2. 计算接触衰减曲线
3. 生成衰减曲线图表
4. 输出统计结果

作者：Claude Code Assistant
日期：2025-08-29
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import warnings
warnings.filterwarnings('ignore')

try:
    import cooler
    import cooltools
except ImportError as e:
    print(f"错误：缺少必要的依赖包: {e}")
    print("请安装：pip install cooler cooltools")
    sys.exit(1)

# 配置matplotlib中文显示
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class ContactDecayAnalyzer:
    """接触衰减曲线分析器"""
    
    def __init__(self, cool_path: str, output_dir: str, log_level: str = 'INFO'):
        """
        初始化分析器
        
        Args:
            cool_path: .cool文件路径
            output_dir: 输出目录路径
            log_level: 日志级别
        """
        self.cool_path = Path(cool_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 初始化属性（必须在_setup_logging之前）
        self.sample_name = self.cool_path.stem
        self.cool_data = None
        self.matrix = None
        self.decay_profile = None
        
        # 设置日志
        self._setup_logging(log_level)
        
    def _setup_logging(self, log_level: str):
        """设置日志配置"""
        log_file = self.output_dir / f"{self.sample_name}_decay_analysis.log"
        logging.basicConfig(
            level=getattr(logging, log_level),
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def load_cool_data(self, resolution: Optional[int] = None) -> bool:
        """
        加载.cool数据
        
        Args:
            resolution: 分辨率，如果为None则使用文件默认分辨率
            
        Returns:
            bool: 是否加载成功
        """
        try:
            self.logger.info(f"开始加载 {self.cool_path}")
            
            if not self.cool_path.exists():
                raise FileNotFoundError(f"文件不存在: {self.cool_path}")
                
            # 加载cooler数据
            self.cool_data = cooler.Cooler(str(self.cool_path))
            
            # 获取基本信息
            self.logger.info(f"分辨率: {self.cool_data.binsize}")
            self.logger.info(f"染色体: {list(self.cool_data.chromnames)}")
            self.logger.info(f"样本: {self.sample_name}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"加载数据失败: {e}")
            return False
    
    def extract_contact_matrix(self, chrom: str = None, start: int = None, 
                              end: int = None, balance: bool = False) -> bool:
        """
        提取接触矩阵
        
        Args:
            chrom: 染色体名称，如果为None则使用第一条染色体
            start: 起始位置
            end: 结束位置  
            balance: 是否使用平衡化数据（默认False，因为imputed数据通常没有balance权重）
            
        Returns:
            bool: 是否提取成功
        """
        try:
            if chrom:
                self.logger.info(f"提取染色体 {chrom} 的接触矩阵")
                self.matrix = self.cool_data.matrix(balance=balance).fetch(chrom)
            else:
                self.logger.info("提取第一条染色体的接触矩阵")
                # 对于全基因组分析，选择第一条主要染色体
                first_chrom = self.cool_data.chromnames[0]
                self.logger.info(f"使用染色体: {first_chrom}")
                self.matrix = self.cool_data.matrix(balance=balance).fetch(first_chrom)
                
            # 替换NaN值为0
            self.matrix = np.nan_to_num(self.matrix, nan=0.0, posinf=0.0, neginf=0.0)
            
            self.logger.info(f"矩阵大小: {self.matrix.shape}")
            self.logger.info(f"非零元素数量: {np.count_nonzero(self.matrix)}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"提取接触矩阵失败: {e}")
            return False
    
    def calculate_decay_profile(self, max_distance: Optional[int] = None, 
                               log_bins: bool = True, n_bins: int = 50) -> bool:
        """
        计算接触衰减曲线
        
        Args:
            max_distance: 最大距离（以bin为单位）
            log_bins: 是否使用对数分箱
            n_bins: 分箱数量
            
        Returns:
            bool: 是否计算成功
        """
        try:
            if self.matrix is None:
                raise ValueError("请先加载接触矩阵")
                
            self.logger.info("开始计算接触衰减曲线")
            
            # 获取矩阵大小
            n = self.matrix.shape[0]
            
            if max_distance is None:
                max_distance = n // 4  # 默认使用矩阵大小的1/4
                
            # 创建距离和接触强度数组
            distances = []
            contacts = []
            
            # 计算不同距离的接触强度
            for d in range(1, min(max_distance, n)):
                # 获取距离为d的所有对角线元素（上三角矩阵）
                if d < n:
                    # 使用numpy的diag函数更高效地获取对角线
                    diag_values = np.diag(self.matrix, k=d)
                    # 过滤有效值
                    valid_mask = ~np.isnan(diag_values) & (diag_values > 0)
                    valid_values = diag_values[valid_mask]
                    
                    if len(valid_values) > 0:
                        distances.append(d * self.cool_data.binsize)  # 转换为物理距离
                        contacts.append(np.mean(valid_values))
                        # 记录详细信息用于调试
                        if d <= 5:  # 只记录前几个距离的详情
                            self.logger.debug(f"距离 {d}: {len(valid_values)} 个有效值, 平均值: {np.mean(valid_values):.6f}")
            
            if not distances:
                raise ValueError("未找到有效的接触数据")
                
            # 创建结果DataFrame
            self.decay_profile = pd.DataFrame({
                'distance': distances,
                'contact_frequency': contacts,
                'distance_kb': np.array(distances) / 1000,  # 转换为kb
                'log_distance': np.log10(distances),
                'log_contact': np.log10(contacts)
            })
            
            self.logger.info(f"成功计算 {len(self.decay_profile)} 个距离点的衰减曲线")
            
            return True
            
        except Exception as e:
            self.logger.error(f"计算衰减曲线失败: {e}")
            return False
    
    def plot_decay_curve(self, save_plots: bool = True, show_plots: bool = False) -> Dict[str, str]:
        """
        绘制接触衰减曲线
        
        Args:
            save_plots: 是否保存图片
            show_plots: 是否显示图片
            
        Returns:
            Dict: 保存的图片路径
        """
        try:
            if self.decay_profile is None:
                raise ValueError("请先计算接触衰减曲线")
                
            plot_paths = {}
            
            # 设置图片样式
            plt.style.use('seaborn-v0_8')
            fig_size = (12, 8)
            
            # 1. 线性坐标图
            fig, ax = plt.subplots(figsize=fig_size)
            ax.plot(self.decay_profile['distance_kb'], self.decay_profile['contact_frequency'], 
                   'o-', markersize=3, linewidth=1.5, color='#1f77b4', alpha=0.8)
            ax.set_xlabel('基因组距离 (kb)', fontsize=12)
            ax.set_ylabel('接触频率', fontsize=12)
            ax.set_title(f'{self.sample_name} - 接触衰减曲线 (线性坐标)', fontsize=14)
            ax.grid(True, alpha=0.3)
            
            if save_plots:
                linear_path = self.output_dir / f"{self.sample_name}_decay_linear.png"
                plt.savefig(linear_path, dpi=300, bbox_inches='tight')
                plot_paths['linear'] = str(linear_path)
                
            if show_plots:
                plt.show()
            else:
                plt.close()
            
            # 2. 双对数坐标图
            fig, ax = plt.subplots(figsize=fig_size)
            ax.loglog(self.decay_profile['distance_kb'], self.decay_profile['contact_frequency'], 
                     'o-', markersize=3, linewidth=1.5, color='#ff7f0e', alpha=0.8)
            ax.set_xlabel('基因组距离 (kb)', fontsize=12)
            ax.set_ylabel('接触频率', fontsize=12)
            ax.set_title(f'{self.sample_name} - 接触衰减曲线 (双对数坐标)', fontsize=14)
            ax.grid(True, alpha=0.3)
            
            # 添加幂律拟合线
            try:
                # 简单的线性拟合在对数空间
                log_dist = np.log10(self.decay_profile['distance_kb'])
                log_contact = np.log10(self.decay_profile['contact_frequency'])
                
                # 去除无穷大值
                valid_mask = np.isfinite(log_dist) & np.isfinite(log_contact)
                if np.sum(valid_mask) > 2:
                    coeffs = np.polyfit(log_dist[valid_mask], log_contact[valid_mask], 1)
                    slope, intercept = coeffs
                    
                    # 绘制拟合线
                    x_fit = np.logspace(np.log10(self.decay_profile['distance_kb'].min()), 
                                       np.log10(self.decay_profile['distance_kb'].max()), 100)
                    y_fit = 10**(slope * np.log10(x_fit) + intercept)
                    ax.plot(x_fit, y_fit, '--', color='red', alpha=0.7, 
                           label=f'幂律拟合: y ∝ x^{slope:.2f}')
                    ax.legend()
                    
                    self.logger.info(f"幂律拟合斜率: {slope:.3f}")
            except:
                self.logger.warning("幂律拟合失败")
            
            if save_plots:
                loglog_path = self.output_dir / f"{self.sample_name}_decay_loglog.png"
                plt.savefig(loglog_path, dpi=300, bbox_inches='tight')
                plot_paths['loglog'] = str(loglog_path)
                
            if show_plots:
                plt.show()
            else:
                plt.close()
            
            # 3. 综合对比图
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            
            # 左图：线性
            ax1.plot(self.decay_profile['distance_kb'], self.decay_profile['contact_frequency'], 
                    'o-', markersize=3, linewidth=1.5, color='#1f77b4', alpha=0.8)
            ax1.set_xlabel('基因组距离 (kb)')
            ax1.set_ylabel('接触频率')
            ax1.set_title('线性坐标')
            ax1.grid(True, alpha=0.3)
            
            # 右图：双对数
            ax2.loglog(self.decay_profile['distance_kb'], self.decay_profile['contact_frequency'], 
                      'o-', markersize=3, linewidth=1.5, color='#ff7f0e', alpha=0.8)
            ax2.set_xlabel('基因组距离 (kb)')
            ax2.set_ylabel('接触频率')
            ax2.set_title('双对数坐标')
            ax2.grid(True, alpha=0.3)
            
            plt.suptitle(f'{self.sample_name} - 接触衰减曲线对比', fontsize=16)
            plt.tight_layout()
            
            if save_plots:
                combined_path = self.output_dir / f"{self.sample_name}_decay_combined.png"
                plt.savefig(combined_path, dpi=300, bbox_inches='tight')
                plot_paths['combined'] = str(combined_path)
                
            if show_plots:
                plt.show()
            else:
                plt.close()
            
            self.logger.info(f"成功生成 {len(plot_paths)} 个图表")
            return plot_paths
            
        except Exception as e:
            self.logger.error(f"绘制图表失败: {e}")
            return {}
    
    def save_results(self) -> str:
        """
        保存分析结果
        
        Returns:
            str: 输出文件路径
        """
        try:
            if self.decay_profile is None:
                raise ValueError("请先计算接触衰减曲线")
                
            # 保存详细结果
            output_file = self.output_dir / f"{self.sample_name}_decay_profile.csv"
            self.decay_profile.to_csv(output_file, index=False)
            
            # 计算摘要统计
            summary = {
                'sample_name': self.sample_name,
                'resolution': self.cool_data.binsize,
                'total_points': len(self.decay_profile),
                'max_distance_kb': self.decay_profile['distance_kb'].max(),
                'min_contact_freq': self.decay_profile['contact_frequency'].min(),
                'max_contact_freq': self.decay_profile['contact_frequency'].max(),
                'mean_contact_freq': self.decay_profile['contact_frequency'].mean()
            }
            
            # 保存摘要
            summary_file = self.output_dir / f"{self.sample_name}_summary.json"
            import json
            with open(summary_file, 'w', encoding='utf-8') as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
            
            self.logger.info(f"结果已保存到: {output_file}")
            self.logger.info(f"摘要已保存到: {summary_file}")
            
            return str(output_file)
            
        except Exception as e:
            self.logger.error(f"保存结果失败: {e}")
            return ""
    
    def run_complete_analysis(self, chrom: str = None, max_distance: int = None, 
                             save_plots: bool = True, show_plots: bool = False) -> Dict:
        """
        运行完整分析流程
        
        Args:
            chrom: 分析的染色体，None为全基因组
            max_distance: 最大分析距离
            save_plots: 是否保存图片
            show_plots: 是否显示图片
            
        Returns:
            Dict: 分析结果摘要
        """
        try:
            self.logger.info(f"开始完整分析流程: {self.sample_name}")
            
            # 1. 加载数据
            if not self.load_cool_data():
                return {'success': False, 'error': '数据加载失败'}
            
            # 2. 提取接触矩阵
            if not self.extract_contact_matrix(chrom=chrom):
                return {'success': False, 'error': '矩阵提取失败'}
            
            # 3. 计算衰减曲线
            if not self.calculate_decay_profile(max_distance=max_distance):
                return {'success': False, 'error': '衰减曲线计算失败'}
            
            # 4. 生成图表
            plot_paths = self.plot_decay_curve(save_plots=save_plots, show_plots=show_plots)
            
            # 5. 保存结果
            result_file = self.save_results()
            
            # 返回结果摘要
            result_summary = {
                'success': True,
                'sample_name': self.sample_name,
                'output_dir': str(self.output_dir),
                'result_file': result_file,
                'plot_paths': plot_paths,
                'total_distance_points': len(self.decay_profile) if self.decay_profile is not None else 0,
                'resolution': self.cool_data.binsize if self.cool_data else None
            }
            
            self.logger.info("完整分析流程成功完成")
            return result_summary
            
        except Exception as e:
            self.logger.error(f"完整分析流程失败: {e}")
            return {'success': False, 'error': str(e)}


def main():
    """主函数，用于命令行调用"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Hi-C接触衰减曲线分析工具')
    parser.add_argument('-i', '--input', required=True, help='输入.cool文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出目录路径')
    parser.add_argument('-c', '--chrom', help='分析指定染色体（默认全基因组）')
    parser.add_argument('-d', '--max-distance', type=int, help='最大分析距离（bin数）')
    parser.add_argument('--show-plots', action='store_true', help='显示图片')
    parser.add_argument('--log-level', default='INFO', help='日志级别')
    
    args = parser.parse_args()
    
    # 创建分析器
    analyzer = ContactDecayAnalyzer(
        cool_path=args.input,
        output_dir=args.output,
        log_level=args.log_level
    )
    
    # 运行分析
    results = analyzer.run_complete_analysis(
        chrom=args.chrom,
        max_distance=args.max_distance,
        save_plots=True,
        show_plots=args.show_plots
    )
    
    # 输出结果
    if results['success']:
        print(f"✅ 分析完成！")
        print(f"📁 输出目录: {results['output_dir']}")
        print(f"📊 数据点数: {results['total_distance_points']}")
        print(f"🔍 分辨率: {results['resolution']}")
    else:
        print(f"❌ 分析失败: {results['error']}")
        sys.exit(1)


if __name__ == "__main__":
    main()