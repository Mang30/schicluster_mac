#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
单个.cool文件接触衰减曲线分析器

功能：
1. 加载.cool文件
2. 计算和分析接触衰减曲线
3. 幂律拟合
4. 生成可视化图表
5. 保存分析结果

作者：Claude Code Assistant
日期：2025-08-29
"""

import json
import logging
from pathlib import Path
import cooler
import cooltools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# --- JSON序列化辅助函数 ---
def convert_numpy_types(obj):
    """
    递归地将字典或列表中的NumPy数值类型转换为Python原生类型。
    这是解决 'Object of type int64 is not JSON serializable' 错误的关键。
    """
    if isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(element) for element in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj

class ContactDecayAnalyzer:
    """单个Hi-C样本接触衰减曲线分析器"""

    def __init__(self, cool_path: str, output_dir: str, log_level: str = 'INFO'):
        self.cool_path = Path(cool_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.sample_name = self.cool_path.stem

        self._setup_logging(log_level)
        
        self.cooler_obj = None
        self.decay_profile = None
        self.fit_params = {}

    def _setup_logging(self, log_level):
        """设置日志"""
        self.logger = logging.getLogger(f"contact_decay_analyzer_{self.sample_name}")
        # 避免重复添加handler
        if not self.logger.handlers:
            self.logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.handlers.append(handler)

    def _load_data(self):
        """加载.cool文件"""
        self.logger.info(f"开始加载 {self.cool_path}")
        self.cooler_obj = cooler.Cooler(str(self.cool_path))
        self.logger.info(f"分辨率: {self.cooler_obj.binsize}")
        self.logger.info(f"染色体: {self.cooler_obj.chromnames}")

    def _calculate_decay_profile(self, chrom: str = None, max_distance: int = None):
        """计算接触衰减曲线"""
        self.logger.info("开始计算接触衰减曲线")
        
        # 如果指定染色体，则创建视图
        if chrom and chrom in self.cooler_obj.chromnames:
            view_df = pd.DataFrame({'name': [chrom]})
            self.logger.info(f"使用染色体: {chrom}")
        else:
            view_df = None # 全基因组
            self.logger.info("使用全基因组")

        expected_df = cooltools.expected.diagsum(
            self.cooler_obj,
            view_df=view_df,
            ignore_diags=2
        )
        
        value_col = 'balanced.avg' if 'balanced.avg.smoothed' in expected_df.columns else 'count.avg'
        
        profile_data = expected_df[['dist', value_col]].dropna()
        profile_data = profile_data[profile_data[value_col] > 0]
        
        # 转换为真实距离
        profile_data['distance_kb'] = profile_data['dist'] * (self.cooler_obj.binsize / 1000)
        profile_data.rename(columns={value_col: 'contact_frequency'}, inplace=True)

        if max_distance:
            profile_data = profile_data[profile_data['dist'] <= max_distance]

        self.decay_profile = profile_data
        self.logger.info(f"成功计算 {len(self.decay_profile)} 个距离点的衰减曲线")

    def _power_law_fit(self):
        """进行幂律拟合"""
        if self.decay_profile is None or self.decay_profile.empty:
            self.logger.warning("衰减曲线数据为空，跳过幂律拟合")
            return

        def power_law(s, a, b):
            return a * (s ** b)

        try:
            # 拟合时忽略非常短的距离
            fit_data = self.decay_profile[self.decay_profile['distance_kb'] > 0]
            params, _ = curve_fit(
                power_law,
                fit_data['distance_kb'],
                fit_data['contact_frequency'],
                maxfev=5000
            )
            self.fit_params = {'a': params[0], 'slope': params[1]}
            self.logger.info(f"幂律拟合斜率: {params[1]:.3f}")
        except Exception as e:
            self.logger.error(f"幂律拟合失败: {e}")
            self.fit_params = {'a': np.nan, 'slope': np.nan}
    
    def _generate_plots(self, save_plots: bool = True):
        """生成可视化图表"""
        if self.decay_profile is None or self.decay_profile.empty:
            return 0
        
        n_plots = 0
        plt.style.use('seaborn-v0_8-whitegrid')

        # 1. 线性坐标图
        plt.figure(figsize=(8, 6))
        plt.plot(self.decay_profile['distance_kb'], self.decay_profile['contact_frequency'], 'o-', markersize=3)
        plt.xlabel('基因组距离 (kb)')
        plt.ylabel('平均接触频率')
        plt.title(f'接触衰减曲线 - {self.sample_name}')
        if save_plots:
            plt.savefig(self.output_dir / "decay_curve_linear.png", dpi=300)
            n_plots += 1
        plt.close()

        # 2. 双对数坐标图及拟合曲线
        plt.figure(figsize=(8, 6))
        plt.loglog(self.decay_profile['distance_kb'], self.decay_profile['contact_frequency'], 'o', markersize=3, label='观测值')
        if 'slope' in self.fit_params and not np.isnan(self.fit_params['slope']):
            fit_y = self.fit_params['a'] * (self.decay_profile['distance_kb'] ** self.fit_params['slope'])
            plt.loglog(self.decay_profile['distance_kb'], fit_y, 'r-', label=f"拟合 (斜率={self.fit_params['slope']:.2f})")
        plt.xlabel('基因组距离 (kb, log scale)')
        plt.ylabel('平均接触频率 (log scale)')
        plt.title(f'接触衰减曲线 (双对数) - {self.sample_name}')
        plt.legend()
        if save_plots:
            plt.savefig(self.output_dir / "decay_curve_loglog.png", dpi=300)
            n_plots += 1
        plt.close()

        self.logger.info(f"成功生成 {n_plots} 个图表")
        return n_plots

    def _save_results(self):
        """保存分析结果到文件"""
        if self.decay_profile is not None:
            self.decay_profile.to_csv(self.output_dir / "decay_profile.csv", index=False)
        
        summary = {
            'sample_name': self.sample_name,
            'cool_file': str(self.cool_path),
            'resolution': self.cooler_obj.binsize if self.cooler_obj else None,
            'total_bins': self.cooler_obj.bins().shape[0] if self.cooler_obj else None,
            'power_law_slope': self.fit_params.get('slope'),
            'data_points': len(self.decay_profile) if self.decay_profile is not None else 0
        }

        # *** 关键修正点：在保存JSON前，对summary字典进行类型转换 ***
        summary = convert_numpy_types(summary)

        try:
            with open(self.output_dir / "analysis_summary.json", 'w') as f:
                json.dump(summary, f, indent=4)
        except TypeError as e:
            self.logger.error(f"保存结果失败: {e}")
            raise # 重新抛出异常，让上层知道
            
        return summary

    def run_complete_analysis(self, chrom: str = None, max_distance: int = None, save_plots: bool = True, show_plots: bool = False):
        """运行完整的分析流程"""
        self.logger.info(f"开始完整分析流程: {self.sample_name}")
        try:
            self._load_data()
            self._calculate_decay_profile(chrom=chrom, max_distance=max_distance)
            self._power_law_fit()
            self._generate_plots(save_plots=save_plots)
            summary = self._save_results()
            summary['success'] = True
            self.logger.info("完整分析流程成功完成")
            return summary
        except Exception as e:
            self.logger.error(f"分析流程失败: {e}")
            return {'success': False, 'error': str(e)}

if __name__ == '__main__':
    # 这个部分用于独立测试单个文件
    parser = argparse.ArgumentParser(description="单个.cool文件接触衰减曲线分析")
    parser.add_argument('-i', '--cool-path', required=True, help=".cool文件路径")
    parser.add_argument('-o', '--output-dir', required=True, help="输出目录")
    parser.add_argument('-c', '--chrom', help="分析指定染色体 (默认全基因组)")
    parser.add_argument('-d', '--max-distance', type=int, help="最大分析距离 (bin数)")
    args = parser.parse_args()

    analyzer = ContactDecayAnalyzer(cool_path=args.cool_path, output_dir=args.output_dir)
    analyzer.run_complete_analysis(chrom=args.chrom, max_distance=args.max_distance)

