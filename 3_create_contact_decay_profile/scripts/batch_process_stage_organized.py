#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量处理stage_organized目录中的接触衰减曲线分析脚本

功能：
1. 自动发现stage_organized目录下的各个stage子目录
2. 批量处理每个stage目录中的.cool文件
3. 生成各个stage的接触衰减曲线
4. 并行处理以提高效率
5. 可以指定处理特定的stage或处理所有stage

作者：Claude Code Assistant
日期：2025-08-31
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# 添加scripts目录到路径
sys.path.insert(0, str(Path(__file__).parent))

try:
    from contact_decay_analyzer import ContactDecayAnalyzer
    print("正确引入ContactDecayAnalyzer模块")
except ImportError:
    print("错误：无法导入ContactDecayAnalyzer，请检查scripts目录或确保其在PYTHONPATH中")
    sys.exit(1)

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# JSON序列化辅助函数
def convert_numpy_types(obj):
    """递归地将字典或列表中的NumPy数值类型转换为Python原生类型"""
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

class StageOrganizedProcessor:
    """Stage组织数据的接触衰减批量处理器"""
    
    def __init__(self, stage_organized_dir: str, output_dir: str, 
                 stages: Optional[List[str]] = None):
        """
        初始化批量处理器
        
        Args:
            stage_organized_dir: stage_organized目录路径
            output_dir: 输出目录路径
            stages: 要处理的发育阶段列表，None表示处理所有发现的阶段
        """
        self.stage_organized_dir = Path(stage_organized_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置日志
        self._setup_logging()
        
        # 发现可用的stage
        self.available_stages = self._discover_stages()
        
        # 确定要处理的stage
        if stages:
            # 验证指定的stage是否存在
            invalid_stages = [s for s in stages if s not in self.available_stages]
            if invalid_stages:
                self.logger.error(f"指定的stage不存在: {invalid_stages}")
                self.logger.info(f"可用的stage: {self.available_stages}")
                sys.exit(1)
            self.stages = stages
        else:
            self.stages = self.available_stages
            
        self.logger.info(f"将要处理的发育阶段: {self.stages}")
        
        # 结果存储
        self.stage_results = {}
        self.all_decay_profiles = {}
        
    def _setup_logging(self):
        """设置日志配置"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = self.output_dir / f"stage_organized_processing_{timestamp}.log"
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def _discover_stages(self) -> List[str]:
        """发现stage_organized目录下的所有stage子目录"""
        stages = []
        if not self.stage_organized_dir.exists():
            self.logger.error(f"Stage organized目录不存在: {self.stage_organized_dir}")
            return stages
            
        for item in self.stage_organized_dir.iterdir():
            if item.is_dir():
                # 检查目录中是否有cool文件
                cool_files = list(item.glob("*.cool"))
                if cool_files:
                    stages.append(item.name)
                    
        stages.sort()
        self.logger.info(f"发现的stage: {stages}")
        return stages
    
    def _get_stage_cool_files(self, stage: str, max_files: Optional[int] = None) -> List[Path]:
        """获取指定stage目录下的所有cool文件"""
        stage_dir = self.stage_organized_dir / stage
        if not stage_dir.exists():
            self.logger.warning(f"Stage目录不存在: {stage_dir}")
            return []
            
        cool_files = list(stage_dir.glob("*.cool"))
        cool_files.sort()
        
        if max_files and len(cool_files) > max_files:
            self.logger.info(f"Stage {stage}: 限制处理文件数从 {len(cool_files)} 到 {max_files}")
            cool_files = cool_files[:max_files]
            
        self.logger.info(f"Stage {stage}: 发现 {len(cool_files)} 个cool文件")
        return cool_files
    
    def process_single_stage(self, stage: str, max_files: Optional[int] = None, 
                           max_workers: int = 4) -> Dict:
        """处理单个stage的所有cool文件"""
        self.logger.info(f"开始处理stage: {stage}")
        
        # 创建stage专用的输出目录
        stage_output_dir = self.output_dir / f"stage_{stage}"
        stage_output_dir.mkdir(parents=True, exist_ok=True)
        
        # 获取cool文件列表
        cool_files = self._get_stage_cool_files(stage, max_files)
        if not cool_files:
            self.logger.warning(f"Stage {stage}: 未找到cool文件")
            return {}
        
        # 初始化结果存储
        stage_results = {
            'stage': stage,
            'total_files': len(cool_files),
            'successful_files': 0,
            'failed_files': 0,
            'cell_decay_profiles': {},
            'summary_statistics': {}
        }
        
        # 并行处理cool文件
        self.logger.info(f"Stage {stage}: 使用 {max_workers} 个worker并行处理")
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_file = {}
            for cool_file in cool_files:
                future = executor.submit(self._process_single_cool_file, cool_file, stage_output_dir)
                future_to_file[future] = cool_file
            
            # 收集结果
            for future in as_completed(future_to_file):
                cool_file = future_to_file[future]
                cell_name = cool_file.stem
                
                try:
                    result = future.result()
                    if result:
                        stage_results['cell_decay_profiles'][cell_name] = result
                        stage_results['successful_files'] += 1
                        if stage_results['successful_files'] % 50 == 0:
                            self.logger.info(f"Stage {stage}: 已完成 {stage_results['successful_files']}/{len(cool_files)} 个文件")
                    else:
                        stage_results['failed_files'] += 1
                        self.logger.warning(f"处理失败: {cell_name}")
                except Exception as exc:
                    stage_results['failed_files'] += 1
                    self.logger.error(f"处理 {cell_name} 时发生异常: {exc}")
        
        # 生成stage汇总统计
        self._generate_stage_summary(stage, stage_results, stage_output_dir)
        
        self.logger.info(f"Stage {stage} 处理完成: 成功 {stage_results['successful_files']}, 失败 {stage_results['failed_files']}")
        return stage_results
    
    def _process_single_cool_file(self, cool_file: Path, output_dir: Path) -> Optional[Dict]:
        """处理单个cool文件"""
        try:
            cell_name = cool_file.stem
            
            # 创建分析器实例
            analyzer = ContactDecayAnalyzer(str(cool_file), str(output_dir))
            
            # 执行分析
            analysis_result = analyzer.run_complete_analysis(save_plots=False, show_plots=False)
            
            if not analysis_result.get('success', False):
                return None
            
            # 从分析器中提取衰减曲线数据
            if analyzer.decay_profile is None or analyzer.decay_profile.empty:
                return None
                
            decay_profile = {
                'distances': analyzer.decay_profile['dist'].tolist(),
                'decay_values': analyzer.decay_profile['contact_frequency'].tolist(),
                'distance_kb': analyzer.decay_profile['distance_kb'].tolist(),
                'total_contacts': int(analyzer.decay_profile['contact_frequency'].sum()),
                'resolution': analyzer.cooler_obj.binsize if analyzer.cooler_obj else None,
                'power_law_slope': analyzer.fit_params.get('slope', np.nan)
            }
            
            # 保存单个细胞的结果
            cell_output_file = output_dir / f"{cell_name}_decay_profile.json"
            with open(cell_output_file, 'w', encoding='utf-8') as f:
                json.dump(convert_numpy_types(decay_profile), f, indent=2, ensure_ascii=False)
            
            # 注意：图表已经在run_complete_analysis中生成，无需重复生成
            
            return decay_profile
            
        except Exception as e:
            logging.error(f"处理 {cool_file} 失败: {e}")
            return None
    
    def _generate_stage_summary(self, stage: str, stage_results: Dict, output_dir: Path):
        """生成stage的汇总统计和图表"""
        try:
            if not stage_results['cell_decay_profiles']:
                self.logger.warning(f"Stage {stage}: 无有效的衰减曲线数据，跳过汇总统计")
                return
            
            # 收集所有细胞的衰减数据
            all_distances = []
            all_decay_values = []
            cell_summaries = []
            
            for cell_name, profile in stage_results['cell_decay_profiles'].items():
                if 'distances' in profile and 'decay_values' in profile:
                    distances = profile['distances']
                    decay_values = profile['decay_values']
                    
                    all_distances.extend(distances)
                    all_decay_values.extend(decay_values)
                    
                    # 收集每个细胞的统计信息
                    cell_summaries.append({
                        'cell_name': cell_name,
                        'mean_decay': np.mean(decay_values),
                        'median_decay': np.median(decay_values),
                        'std_decay': np.std(decay_values),
                        'total_contacts': profile.get('total_contacts', 0)
                    })
            
            # 计算stage级别的统计
            stage_stats = {
                'total_cells': len(cell_summaries),
                'mean_decay': np.mean(all_decay_values) if all_decay_values else 0,
                'median_decay': np.median(all_decay_values) if all_decay_values else 0,
                'std_decay': np.std(all_decay_values) if all_decay_values else 0,
                'total_contacts': sum(cell['total_contacts'] for cell in cell_summaries)
            }
            
            stage_results['summary_statistics'] = stage_stats
            
            # 保存汇总结果
            summary_file = output_dir / f"stage_{stage}_summary.json"
            with open(summary_file, 'w', encoding='utf-8') as f:
                json.dump(convert_numpy_types(stage_results), f, indent=2, ensure_ascii=False)
            
            # 保存细胞级别的汇总表格
            if cell_summaries:
                df = pd.DataFrame(cell_summaries)
                df.to_csv(output_dir / f"stage_{stage}_cell_summary.csv", index=False)
            
            # 生成汇总图表
            self._plot_stage_summary(stage, stage_results, output_dir)
            
        except Exception as e:
            self.logger.error(f"生成Stage {stage} 汇总统计失败: {e}")
    
    def _plot_stage_summary(self, stage: str, stage_results: Dict, output_dir: Path):
        """生成stage的汇总图表"""
        try:
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle(f'Stage {stage} - 接触衰减曲线汇总分析', fontsize=16)
            
            # 收集数据
            decay_means = []
            total_contacts = []
            cell_names = []
            
            for cell_name, profile in stage_results['cell_decay_profiles'].items():
                if 'decay_values' in profile:
                    decay_means.append(np.mean(profile['decay_values']))
                    total_contacts.append(profile.get('total_contacts', 0))
                    cell_names.append(cell_name)
            
            if not decay_means:
                self.logger.warning(f"Stage {stage}: 无数据可绘制")
                return
            
            # 1. 衰减值分布直方图
            axes[0, 0].hist(decay_means, bins=30, alpha=0.7, color='skyblue')
            axes[0, 0].set_title('细胞间平均衰减值分布')
            axes[0, 0].set_xlabel('平均衰减值')
            axes[0, 0].set_ylabel('细胞数量')
            
            # 2. 接触数分布
            axes[0, 1].hist(total_contacts, bins=30, alpha=0.7, color='lightgreen')
            axes[0, 1].set_title('细胞间总接触数分布')
            axes[0, 1].set_xlabel('总接触数')
            axes[0, 1].set_ylabel('细胞数量')
            
            # 3. 衰减值vs接触数散点图
            axes[1, 0].scatter(total_contacts, decay_means, alpha=0.6, color='coral')
            axes[1, 0].set_title('平均衰减值 vs 总接触数')
            axes[1, 0].set_xlabel('总接触数')
            axes[1, 0].set_ylabel('平均衰减值')
            
            # 4. 统计信息文本
            stats_text = f"""Stage {stage} 统计信息:
细胞总数: {stage_results['successful_files']}
平均衰减值: {stage_results['summary_statistics']['mean_decay']:.6f}
衰减值标准差: {stage_results['summary_statistics']['std_decay']:.6f}
总接触数: {stage_results['summary_statistics']['total_contacts']:,}
"""
            axes[1, 1].text(0.1, 0.5, stats_text, fontsize=12, verticalalignment='center')
            axes[1, 1].set_xlim(0, 1)
            axes[1, 1].set_ylim(0, 1)
            axes[1, 1].axis('off')
            
            plt.tight_layout()
            plt.savefig(output_dir / f"stage_{stage}_summary_plot.png", dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            self.logger.error(f"绘制Stage {stage} 汇总图表失败: {e}")
    
    def process_all_stages(self, max_files: Optional[int] = None, max_workers: int = 4):
        """处理所有指定的stage"""
        self.logger.info(f"开始批量处理 {len(self.stages)} 个stage")
        
        for i, stage in enumerate(self.stages, 1):
            self.logger.info(f"处理进度 [{i}/{len(self.stages)}]: Stage {stage}")
            
            stage_results = self.process_single_stage(stage, max_files, max_workers)
            self.stage_results[stage] = stage_results
        
        # 生成跨stage的比较分析
        self._generate_cross_stage_analysis()
        
        self.logger.info("所有stage处理完成！")
    
    def _generate_cross_stage_analysis(self):
        """生成跨stage的比较分析"""
        try:
            # 收集各stage的汇总统计
            comparison_data = []
            for stage, results in self.stage_results.items():
                if 'summary_statistics' in results:
                    stats = results['summary_statistics']
                    comparison_data.append({
                        'stage': stage,
                        'total_cells': stats['total_cells'],
                        'mean_decay': stats['mean_decay'],
                        'median_decay': stats['median_decay'],
                        'std_decay': stats['std_decay'],
                        'total_contacts': stats['total_contacts']
                    })
            
            if not comparison_data:
                self.logger.warning("无数据进行跨stage比较分析")
                return
            
            # 保存比较数据
            df_comparison = pd.DataFrame(comparison_data)
            df_comparison.to_csv(self.output_dir / "cross_stage_comparison.csv", index=False)
            
            # 生成比较图表
            self._plot_cross_stage_comparison(df_comparison)
            
            # 保存完整的结果
            final_results = {
                'processing_summary': {
                    'total_stages': len(self.stages),
                    'processed_stages': list(self.stage_results.keys()),
                    'timestamp': datetime.now().isoformat()
                },
                'stage_results': self.stage_results,
                'cross_stage_comparison': comparison_data
            }
            
            with open(self.output_dir / "complete_analysis_results.json", 'w', encoding='utf-8') as f:
                json.dump(convert_numpy_types(final_results), f, indent=2, ensure_ascii=False)
                
        except Exception as e:
            self.logger.error(f"生成跨stage比较分析失败: {e}")
    
    def _plot_cross_stage_comparison(self, df: pd.DataFrame):
        """绘制跨stage的比较图表"""
        try:
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            fig.suptitle('跨Stage接触衰减曲线比较分析', fontsize=16)
            
            stages = df['stage'].tolist()
            
            # 1. 各stage平均衰减值比较
            axes[0, 0].bar(stages, df['mean_decay'], color='skyblue', alpha=0.7)
            axes[0, 0].set_title('各Stage平均衰减值比较')
            axes[0, 0].set_ylabel('平均衰减值')
            axes[0, 0].tick_params(axis='x', rotation=45)
            
            # 2. 各stage细胞数量比较
            axes[0, 1].bar(stages, df['total_cells'], color='lightgreen', alpha=0.7)
            axes[0, 1].set_title('各Stage细胞数量比较')
            axes[0, 1].set_ylabel('细胞数量')
            axes[0, 1].tick_params(axis='x', rotation=45)
            
            # 3. 各stage总接触数比较
            axes[1, 0].bar(stages, df['total_contacts'], color='coral', alpha=0.7)
            axes[1, 0].set_title('各Stage总接触数比较')
            axes[1, 0].set_ylabel('总接触数')
            axes[1, 0].tick_params(axis='x', rotation=45)
            
            # 4. 衰减值标准差比较
            axes[1, 1].bar(stages, df['std_decay'], color='gold', alpha=0.7)
            axes[1, 1].set_title('各Stage衰减值标准差比较')
            axes[1, 1].set_ylabel('标准差')
            axes[1, 1].tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / "cross_stage_comparison_plot.png", dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            self.logger.error(f"绘制跨stage比较图表失败: {e}")

def main():
    parser = argparse.ArgumentParser(description='批量处理stage_organized目录中的接触衰减曲线分析')
    parser.add_argument('--input', '-i', required=True,
                       help='stage_organized目录路径')
    parser.add_argument('--output', '-o', required=True,
                       help='输出目录路径')
    parser.add_argument('--stages', '-s',
                       help='要处理的stage列表，用逗号分隔，例如：E75,E85,EX05')
    parser.add_argument('--max-workers', '-w', type=int, default=4,
                       help='并行处理的最大worker数量')
    parser.add_argument('--max-files', '-n', type=int,
                       help='每个stage最大处理文件数量（用于测试）')
    
    args = parser.parse_args()
    
    # 解析stages参数
    stages = None
    if args.stages:
        stages = [s.strip() for s in args.stages.split(',')]
    
    # 创建处理器
    processor = StageOrganizedProcessor(
        stage_organized_dir=args.input,
        output_dir=args.output,
        stages=stages
    )
    
    # 执行处理
    processor.process_all_stages(
        max_files=args.max_files,
        max_workers=args.max_workers
    )

if __name__ == "__main__":
    main()