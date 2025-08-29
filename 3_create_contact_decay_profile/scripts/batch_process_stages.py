#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量处理各个发育阶段的接触衰减曲线分析脚本

功能：
1. 自动发现各个stage目录下的.cool文件
2. 批量运行接触衰减曲线分析
3. 生成对比分析结果
4. 并行处理以提高效率

作者：Claude Code Assistant
日期：2025-08-29
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
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# 添加src目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

try:
    from contact_decay_analyzer import ContactDecayAnalyzer
except ImportError:
    print("错误：无法导入ContactDecayAnalyzer，请检查src目录")
    sys.exit(1)

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class StageWiseDecayProcessor:
    """发育阶段接触衰减批量处理器"""
    
    def __init__(self, data_root_dir: str, output_dir: str, 
                 stages: Optional[List[str]] = None):
        """
        初始化批量处理器
        
        Args:
            data_root_dir: 数据根目录路径
            output_dir: 输出目录路径
            stages: 要处理的发育阶段列表，None表示处理所有发现的阶段
        """
        self.data_root_dir = Path(data_root_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置日志
        self._setup_logging()
        
        # 发育阶段
        self.stages = stages if stages else self._discover_stages()
        self.logger.info(f"发现的发育阶段: {self.stages}")
        
        # 结果存储
        self.stage_results = {}
        self.all_decay_profiles = {}
        
    def _setup_logging(self):
        """设置日志配置"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = self.output_dir / f"batch_processing_{timestamp}.log"
        
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
        """自动发现数据目录中的发育阶段"""
        stages = []
        
        try:
            # 查找以E或EX开头的目录
            for item in self.data_root_dir.iterdir():
                if item.is_dir() and (item.name.startswith('E') or item.name.startswith('EX')):
                    # 检查是否包含插补数据
                    impute_dir = item / 'impute' / '100K' / 'chunk0'
                    if impute_dir.exists() and any(impute_dir.glob('*.cool')):
                        stages.append(item.name)
                        
            # 按名称排序
            stages.sort(key=lambda x: (len(x), x))
            return stages
            
        except Exception as e:
            self.logger.error(f"发现发育阶段失败: {e}")
            return []
    
    def _get_stage_cool_files(self, stage: str, max_files: Optional[int] = None) -> List[Path]:
        """
        获取指定阶段的.cool文件列表
        
        Args:
            stage: 发育阶段名称
            max_files: 最大文件数量限制
            
        Returns:
            List[Path]: .cool文件路径列表
        """
        cool_files = []
        
        try:
            stage_dir = self.data_root_dir / stage / 'impute' / '100K' / 'chunk0'
            if not stage_dir.exists():
                self.logger.warning(f"阶段 {stage} 的数据目录不存在: {stage_dir}")
                return []
            
            # 查找所有.cool文件
            cool_files = list(stage_dir.glob('*.cool'))
            
            # 限制文件数量
            if max_files and len(cool_files) > max_files:
                cool_files = cool_files[:max_files]
                self.logger.info(f"阶段 {stage}: 限制处理文件数量为 {max_files}")
            
            self.logger.info(f"阶段 {stage}: 发现 {len(cool_files)} 个.cool文件")
            return cool_files
            
        except Exception as e:
            self.logger.error(f"获取阶段 {stage} 的文件列表失败: {e}")
            return []
    
    def process_single_file(self, cool_file: Path, stage: str, 
                           output_subdir: str, **kwargs) -> Dict:
        """
        处理单个.cool文件
        
        Args:
            cool_file: .cool文件路径
            stage: 发育阶段
            output_subdir: 输出子目录
            **kwargs: 传递给分析器的参数
            
        Returns:
            Dict: 分析结果
        """
        try:
            # 创建输出目录
            file_output_dir = self.output_dir / output_subdir / cool_file.stem
            
            # 创建分析器
            analyzer = ContactDecayAnalyzer(
                cool_path=str(cool_file),
                output_dir=str(file_output_dir),
                log_level='WARNING'  # 减少并行处理时的日志输出
            )
            
            # 运行分析
            results = analyzer.run_complete_analysis(
                chrom=kwargs.get('chrom'),
                max_distance=kwargs.get('max_distance'),
                save_plots=kwargs.get('save_plots', True),
                show_plots=False
            )
            
            # 添加阶段信息
            results['stage'] = stage
            results['cool_file'] = str(cool_file)
            
            return results
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'stage': stage,
                'cool_file': str(cool_file)
            }
    
    def process_stage(self, stage: str, max_files: Optional[int] = None, 
                     parallel: bool = True, max_workers: int = 4, **kwargs) -> Dict:
        """
        处理单个发育阶段的所有.cool文件
        
        Args:
            stage: 发育阶段名称
            max_files: 最大处理文件数量
            parallel: 是否使用并行处理
            max_workers: 并行处理的最大worker数量
            **kwargs: 传递给分析的参数
            
        Returns:
            Dict: 阶段处理结果摘要
        """
        self.logger.info(f"开始处理发育阶段: {stage}")
        
        # 获取文件列表
        cool_files = self._get_stage_cool_files(stage, max_files)
        if not cool_files:
            return {'stage': stage, 'success': False, 'error': '未找到.cool文件'}
        
        # 创建阶段输出目录
        stage_output_dir = f"stage_{stage}"
        
        # 处理文件
        results = []
        successful_count = 0
        failed_count = 0
        
        if parallel and len(cool_files) > 1:
            # 并行处理
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # 提交任务
                future_to_file = {
                    executor.submit(
                        self.process_single_file, 
                        cool_file, stage, stage_output_dir, **kwargs
                    ): cool_file 
                    for cool_file in cool_files
                }
                
                # 收集结果
                for future in as_completed(future_to_file):
                    cool_file = future_to_file[future]
                    try:
                        result = future.result()
                        results.append(result)
                        
                        if result.get('success', False):
                            successful_count += 1
                            self.logger.info(f"✅ {cool_file.stem}: 处理成功")
                        else:
                            failed_count += 1
                            self.logger.error(f"❌ {cool_file.stem}: {result.get('error', '未知错误')}")
                            
                    except Exception as e:
                        failed_count += 1
                        self.logger.error(f"❌ {cool_file.stem}: 执行异常 - {e}")
                        results.append({
                            'success': False,
                            'error': str(e),
                            'stage': stage,
                            'cool_file': str(cool_file)
                        })
        else:
            # 串行处理
            for cool_file in cool_files:
                result = self.process_single_file(cool_file, stage, stage_output_dir, **kwargs)
                results.append(result)
                
                if result.get('success', False):
                    successful_count += 1
                    self.logger.info(f"✅ {cool_file.stem}: 处理成功")
                else:
                    failed_count += 1
                    self.logger.error(f"❌ {cool_file.stem}: {result.get('error', '未知错误')}")
        
        # 汇总结果
        stage_summary = {
            'stage': stage,
            'total_files': len(cool_files),
            'successful': successful_count,
            'failed': failed_count,
            'success_rate': successful_count / len(cool_files) if cool_files else 0,
            'results': results,
            'output_dir': str(self.output_dir / stage_output_dir)
        }
        
        # 保存阶段摘要
        summary_file = self.output_dir / f"stage_{stage}_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(stage_summary, f, indent=2, ensure_ascii=False, default=str)
        
        self.logger.info(f"阶段 {stage} 处理完成: {successful_count}/{len(cool_files)} 成功")
        return stage_summary
    
    def collect_decay_profiles(self) -> pd.DataFrame:
        """收集所有成功分析的衰减曲线数据"""
        all_profiles = []
        
        for stage in self.stages:
            stage_dir = self.output_dir / f"stage_{stage}"
            if not stage_dir.exists():
                continue
                
            # 查找所有衰减曲线文件
            for profile_file in stage_dir.glob('*/*/decay_profile.csv'):
                try:
                    df = pd.read_csv(profile_file)
                    df['stage'] = stage
                    df['sample'] = profile_file.parent.name
                    all_profiles.append(df)
                except Exception as e:
                    self.logger.warning(f"无法读取文件 {profile_file}: {e}")
        
        if all_profiles:
            combined_df = pd.concat(all_profiles, ignore_index=True)
            
            # 保存合并的数据
            output_file = self.output_dir / 'all_stages_decay_profiles.csv'
            combined_df.to_csv(output_file, index=False)
            self.logger.info(f"合并的衰减曲线数据已保存到: {output_file}")
            
            return combined_df
        else:
            self.logger.warning("未找到任何衰减曲线数据")
            return pd.DataFrame()
    
    def create_comparative_plots(self, combined_df: pd.DataFrame) -> List[str]:
        """创建各阶段间的对比图表"""
        if combined_df.empty:
            return []
        
        plot_files = []
        
        try:
            # 设置颜色方案
            stage_colors = plt.cm.Set3(np.linspace(0, 1, len(self.stages)))
            color_map = dict(zip(self.stages, stage_colors))
            
            # 1. 各阶段平均衰减曲线对比
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            
            # 计算每个阶段的平均衰减曲线
            for stage in self.stages:
                stage_data = combined_df[combined_df['stage'] == stage]
                if stage_data.empty:
                    continue
                    
                # 按距离分组计算平均值
                avg_decay = stage_data.groupby('distance_kb')['contact_frequency'].mean()
                
                if not avg_decay.empty:
                    # 线性图
                    ax1.plot(avg_decay.index, avg_decay.values, 
                            'o-', markersize=3, linewidth=2, 
                            label=f'阶段 {stage}', color=color_map[stage], alpha=0.8)
                    
                    # 双对数图
                    ax2.loglog(avg_decay.index, avg_decay.values, 
                              'o-', markersize=3, linewidth=2, 
                              label=f'阶段 {stage}', color=color_map[stage], alpha=0.8)
            
            # 设置图表
            ax1.set_xlabel('基因组距离 (kb)')
            ax1.set_ylabel('平均接触频率')
            ax1.set_title('各发育阶段接触衰减对比 (线性坐标)')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            ax2.set_xlabel('基因组距离 (kb)')
            ax2.set_ylabel('平均接触频率')
            ax2.set_title('各发育阶段接触衰减对比 (双对数坐标)')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # 保存图片
            plot_file = self.output_dir / 'stage_comparison_decay_curves.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plot_files.append(str(plot_file))
            plt.close()
            
            # 2. 热图展示不同距离下各阶段的接触强度
            if len(self.stages) > 1:
                # 准备热图数据
                pivot_data = []
                distance_points = sorted(combined_df['distance_kb'].unique())[:20]  # 取前20个距离点
                
                for stage in self.stages:
                    stage_data = combined_df[combined_df['stage'] == stage]
                    stage_avg = []
                    
                    for dist in distance_points:
                        dist_data = stage_data[abs(stage_data['distance_kb'] - dist) < 10]  # 允许一定误差
                        if not dist_data.empty:
                            stage_avg.append(dist_data['contact_frequency'].mean())
                        else:
                            stage_avg.append(np.nan)
                    
                    pivot_data.append(stage_avg)
                
                # 创建热图
                fig, ax = plt.subplots(figsize=(12, 6))
                heatmap_data = np.array(pivot_data)
                
                # 使用对数变换
                heatmap_data_log = np.log10(heatmap_data + 1e-10)  # 避免log(0)
                
                im = ax.imshow(heatmap_data_log, aspect='auto', cmap='viridis')
                
                # 设置标签
                ax.set_xticks(range(len(distance_points)))
                ax.set_xticklabels([f'{int(d)}kb' for d in distance_points], rotation=45)
                ax.set_yticks(range(len(self.stages)))
                ax.set_yticklabels(self.stages)
                
                ax.set_xlabel('基因组距离')
                ax.set_ylabel('发育阶段')
                ax.set_title('各发育阶段接触强度热图 (log10变换)')
                
                # 添加颜色条
                plt.colorbar(im, ax=ax, label='log10(接触频率)')
                plt.tight_layout()
                
                # 保存图片
                heatmap_file = self.output_dir / 'stage_contact_heatmap.png'
                plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
                plot_files.append(str(heatmap_file))
                plt.close()
            
            self.logger.info(f"生成对比图表: {len(plot_files)} 个")
            return plot_files
            
        except Exception as e:
            self.logger.error(f"生成对比图表失败: {e}")
            return []
    
    def run_batch_processing(self, max_files_per_stage: Optional[int] = None,
                           parallel: bool = True, max_workers: int = 4,
                           **analysis_kwargs) -> Dict:
        """
        运行批量处理流程
        
        Args:
            max_files_per_stage: 每个阶段最大处理文件数
            parallel: 是否使用并行处理
            max_workers: 最大worker数量
            **analysis_kwargs: 传递给分析的参数
            
        Returns:
            Dict: 批量处理结果摘要
        """
        start_time = datetime.now()
        self.logger.info("开始批量处理所有发育阶段")
        
        # 处理每个阶段
        all_results = {}
        total_successful = 0
        total_files = 0
        
        for stage in self.stages:
            stage_result = self.process_stage(
                stage=stage,
                max_files=max_files_per_stage,
                parallel=parallel,
                max_workers=max_workers,
                **analysis_kwargs
            )
            
            all_results[stage] = stage_result
            total_successful += stage_result.get('successful', 0)
            total_files += stage_result.get('total_files', 0)
        
        # 收集衰减曲线数据
        combined_df = self.collect_decay_profiles()
        
        # 生成对比图表
        comparison_plots = []
        if not combined_df.empty:
            comparison_plots = self.create_comparative_plots(combined_df)
        
        # 生成最终报告
        end_time = datetime.now()
        processing_time = (end_time - start_time).total_seconds()
        
        final_summary = {
            'start_time': start_time.isoformat(),
            'end_time': end_time.isoformat(),
            'processing_time_seconds': processing_time,
            'stages_processed': len(self.stages),
            'total_files': total_files,
            'total_successful': total_successful,
            'overall_success_rate': total_successful / total_files if total_files > 0 else 0,
            'stage_results': all_results,
            'output_directory': str(self.output_dir),
            'comparison_plots': comparison_plots,
            'combined_data_points': len(combined_df) if not combined_df.empty else 0
        }
        
        # 保存最终报告
        report_file = self.output_dir / 'batch_processing_report.json'
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(final_summary, f, indent=2, ensure_ascii=False, default=str)
        
        self.logger.info(f"批量处理完成！")
        self.logger.info(f"总耗时: {processing_time:.1f} 秒")
        self.logger.info(f"处理文件: {total_successful}/{total_files} 成功")
        self.logger.info(f"最终报告: {report_file}")
        
        return final_summary


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='批量处理发育阶段接触衰减曲线分析')
    
    # 必需参数
    parser.add_argument('-i', '--input', required=True, 
                       help='输入数据根目录 (包含各个stage子目录)')
    parser.add_argument('-o', '--output', required=True, 
                       help='输出目录')
    
    # 可选参数
    parser.add_argument('-s', '--stages', nargs='+', 
                       help='指定要处理的发育阶段 (默认处理所有发现的阶段)')
    parser.add_argument('-n', '--max-files', type=int, 
                       help='每个阶段最大处理文件数量')
    parser.add_argument('-w', '--max-workers', type=int, default=4,
                       help='并行处理的最大worker数量')
    parser.add_argument('--no-parallel', action='store_true',
                       help='禁用并行处理')
    
    # 分析参数
    parser.add_argument('-c', '--chrom', help='分析指定染色体')
    parser.add_argument('-d', '--max-distance', type=int, help='最大分析距离')
    parser.add_argument('--no-plots', action='store_true', help='不生成图表')
    
    args = parser.parse_args()
    
    # 创建处理器
    processor = StageWiseDecayProcessor(
        data_root_dir=args.input,
        output_dir=args.output,
        stages=args.stages
    )
    
    # 运行批量处理
    results = processor.run_batch_processing(
        max_files_per_stage=args.max_files,
        parallel=not args.no_parallel,
        max_workers=args.max_workers,
        chrom=args.chrom,
        max_distance=args.max_distance,
        save_plots=not args.no_plots
    )
    
    # 输出结果摘要
    print(f"\n🎉 批量处理完成！")
    print(f"📊 处理了 {results['stages_processed']} 个发育阶段")
    print(f"📁 成功分析 {results['total_successful']}/{results['total_files']} 个文件")
    print(f"⏱️  总耗时: {results['processing_time_seconds']:.1f} 秒")
    print(f"📈 成功率: {results['overall_success_rate']*100:.1f}%")
    print(f"📂 输出目录: {results['output_directory']}")
    
    if results['comparison_plots']:
        print(f"📊 生成对比图表: {len(results['comparison_plots'])} 个")


if __name__ == "__main__":
    main()