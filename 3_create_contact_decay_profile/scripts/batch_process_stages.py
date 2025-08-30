#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量处理各个发育阶段的接触衰减曲线分析脚本 - 优化版

功能：
1. 自动发现各个stage和resolution目录下的.cool文件
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
import numpy as np # 引入numpy用于类型检查
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# 添加src目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

try:
    from contact_decay_analyzer import ContactDecayAnalyzer
    print("正确引入ContactDecayAnalyzer模块")
except ImportError:
    print("错误：无法导入ContactDecayAnalyzer，请检查src目录或确保其在PYTHONPATH中")
    sys.exit(1)

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# --- 新增：JSON序列化辅助函数 ---
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

class StageWiseDecayProcessor:
    """发育阶段接触衰减批量处理器"""
    
    def __init__(self, data_root_dir: str, output_dir: str, 
                 resolution: str, stages: Optional[List[str]] = None):
        """
        初始化批量处理器
        
        Args:
            data_root_dir: 数据根目录路径
            output_dir: 输出目录路径
            resolution: 要处理的分辨率目录名, 例如 '100K'
            stages: 要处理的发育阶段列表，None表示处理所有发现的阶段
        """
        self.data_root_dir = Path(data_root_dir)
        self.output_dir = Path(output_dir)
        self.resolution = resolution
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 设置日志
        self._setup_logging()
        
        # 发育阶段
        self.stages = stages if stages else self._discover_stages()
        self.logger.info(f"将要处理的发育阶段: {self.stages}")
        
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
                    # 检查是否包含指定分辨率的插补数据
                    impute_dir = item / 'impute' / self.resolution / 'chunk0'
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
            # 使用 self.resolution 来构建路径
            stage_dir = self.data_root_dir / stage / 'impute' / self.resolution
            # 搜索所有 chunk* 目录下的 .cool 文件
            cool_files = sorted(stage_dir.glob('chunk*/*.cool'))
            
            if not cool_files:
                self.logger.warning(f"阶段 {stage} 在 {self.resolution} 分辨率下未找到.cool文件: {stage_dir}")
                return []
            
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
        """
        try:
            file_output_dir = self.output_dir / output_subdir / cool_file.stem
            analyzer = ContactDecayAnalyzer(
                cool_path=str(cool_file),
                output_dir=str(file_output_dir),
                log_level='WARNING'
            )
            results = analyzer.run_complete_analysis(
                chrom=kwargs.get('chrom'),
                max_distance=kwargs.get('max_distance'),
                save_plots=kwargs.get('save_plots', True),
                show_plots=False
            )
            results['stage'] = stage
            results['cool_file'] = str(cool_file)
            return results
        except Exception as e:
            return {'success': False, 'error': str(e), 'stage': stage, 'cool_file': str(cool_file)}
    
    def process_stage(self, stage: str, max_files: Optional[int] = None, 
                     parallel: bool = True, max_workers: int = 4, **kwargs) -> Dict:
        """
        处理单个发育阶段的所有.cool文件
        """
        self.logger.info(f"开始处理发育阶段: {stage} @ {self.resolution}")
        cool_files = self._get_stage_cool_files(stage, max_files)
        if not cool_files:
            return {'stage': stage, 'success': False, 'error': '未找到.cool文件'}
        
        stage_output_dir = f"stage_{stage}_{self.resolution}"
        results, successful_count, failed_count = [], 0, 0
        
        if parallel and len(cool_files) > 1:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                future_to_file = {executor.submit(self.process_single_file, cool_file, stage, stage_output_dir, **kwargs): cool_file for cool_file in cool_files}
                for future in as_completed(future_to_file):
                    cool_file = future_to_file[future]
                    try:
                        result = future.result()
                        results.append(result)
                        if result.get('success', False): successful_count += 1
                        else: failed_count += 1; self.logger.error(f"❌ {cool_file.stem}: {result.get('error', '未知错误')}")
                    except Exception as e:
                        failed_count += 1; self.logger.error(f"❌ {cool_file.stem}: 执行异常 - {e}")
                        results.append({'success': False, 'error': str(e), 'stage': stage, 'cool_file': str(cool_file)})
        else:
            for cool_file in cool_files:
                result = self.process_single_file(cool_file, stage, stage_output_dir, **kwargs)
                results.append(result)
                if result.get('success', False): successful_count += 1
                else: failed_count += 1; self.logger.error(f"❌ {cool_file.stem}: {result.get('error', '未知错误')}")

        stage_summary = {
            'stage': stage, 'resolution': self.resolution, 'total_files': len(cool_files),
            'successful': successful_count, 'failed': failed_count,
            'success_rate': successful_count / len(cool_files) if cool_files else 0,
            'results': results, 'output_dir': str(self.output_dir / stage_output_dir)
        }
        
        # *** 修正点：在保存JSON前进行类型转换 ***
        stage_summary = convert_numpy_types(stage_summary)
        
        summary_file = self.output_dir / f"stage_{stage}_{self.resolution}_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f: json.dump(stage_summary, f, indent=2, ensure_ascii=False, default=str)
        self.logger.info(f"阶段 {stage} 处理完成: {successful_count}/{len(cool_files)} 成功")
        return stage_summary
    
    # ... (其他函数 collect_decay_profiles, create_comparative_plots, run_batch_processing 保持不变) ...
    def collect_decay_profiles(self) -> pd.DataFrame:
        """收集所有成功分析的衰减曲线数据"""
        all_profiles = []
        for stage in self.stages:
            stage_dir = self.output_dir / f"stage_{stage}_{self.resolution}"
            if not stage_dir.exists(): continue
            for profile_file in stage_dir.glob('*/*/decay_profile.csv'):
                try:
                    df = pd.read_csv(profile_file)
                    df['stage'] = stage
                    df['resolution'] = self.resolution
                    df['sample'] = profile_file.parent.name
                    all_profiles.append(df)
                except Exception as e: self.logger.warning(f"无法读取文件 {profile_file}: {e}")
        if all_profiles:
            combined_df = pd.concat(all_profiles, ignore_index=True)
            output_file = self.output_dir / f'all_stages_{self.resolution}_decay_profiles.csv'
            combined_df.to_csv(output_file, index=False)
            self.logger.info(f"合并的衰减曲线数据已保存到: {output_file}")
            return combined_df
        else:
            self.logger.warning("未找到任何衰减曲线数据")
            return pd.DataFrame()

    def create_comparative_plots(self, combined_df: pd.DataFrame) -> List[str]:
        """创建各阶段间的对比图表"""
        # (此函数无需修改，因为它依赖于DataFrame的内容)
        if combined_df.empty: return []
        plot_files = []
        try:
            stage_colors = plt.cm.Set3(np.linspace(0, 1, len(self.stages)))
            color_map = dict(zip(self.stages, stage_colors))
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
            for stage in self.stages:
                stage_data = combined_df[combined_df['stage'] == stage]
                if stage_data.empty: continue
                avg_decay = stage_data.groupby('distance_kb')['contact_frequency'].mean()
                if not avg_decay.empty:
                    ax1.plot(avg_decay.index, avg_decay.values, 'o-', markersize=3, linewidth=2, label=f'阶段 {stage}', color=color_map[stage], alpha=0.8)
                    ax2.loglog(avg_decay.index, avg_decay.values, 'o-', markersize=3, linewidth=2, label=f'阶段 {stage}', color=color_map[stage], alpha=0.8)
            ax1.set_xlabel('基因组距离 (kb)'); ax1.set_ylabel('平均接触频率'); ax1.set_title('各发育阶段接触衰减对比 (线性坐标)'); ax1.legend(); ax1.grid(True, alpha=0.3)
            ax2.set_xlabel('基因组距离 (kb)'); ax2.set_ylabel('平均接触频率'); ax2.set_title('各发育阶段接触衰减对比 (双对数坐标)'); ax2.legend(); ax2.grid(True, alpha=0.3)
            plt.tight_layout()
            plot_file = self.output_dir / f'stage_comparison_decay_curves_{self.resolution}.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight'); plot_files.append(str(plot_file)); plt.close()
            self.logger.info(f"生成对比图表: {len(plot_files)} 个")
            return plot_files
        except Exception as e:
            self.logger.error(f"生成对比图表失败: {e}"); return []

    def run_batch_processing(self, max_files_per_stage: Optional[int] = None, parallel: bool = True, max_workers: int = 4, **analysis_kwargs) -> Dict:
        """运行批量处理流程"""
        # (此函数无需修改)
        start_time = datetime.now()
        self.logger.info("开始批量处理所有发育阶段")
        all_results, total_successful, total_files = {}, 0, 0
        for stage in self.stages:
            stage_result = self.process_stage(stage=stage, max_files=max_files_per_stage, parallel=parallel, max_workers=max_workers, **analysis_kwargs)
            all_results[stage] = stage_result
            total_successful += stage_result.get('successful', 0); total_files += stage_result.get('total_files', 0)
        combined_df = self.collect_decay_profiles()
        comparison_plots = []
        if not combined_df.empty: comparison_plots = self.create_comparative_plots(combined_df)
        end_time = datetime.now(); processing_time = (end_time - start_time).total_seconds()
        final_summary = {
            'start_time': start_time.isoformat(), 'end_time': end_time.isoformat(), 'processing_time_seconds': processing_time,
            'stages_processed': len(self.stages), 'total_files': total_files, 'total_successful': total_successful,
            'overall_success_rate': total_successful / total_files if total_files > 0 else 0, 'stage_results': all_results,
            'output_directory': str(self.output_dir), 'comparison_plots': comparison_plots,
            'combined_data_points': len(combined_df) if not combined_df.empty else 0
        }
        
        # *** 修正点：在保存JSON前进行类型转换 ***
        final_summary = convert_numpy_types(final_summary)

        report_file = self.output_dir / f'batch_processing_report_{self.resolution}.json'
        with open(report_file, 'w', encoding='utf-8') as f: json.dump(final_summary, f, indent=2, ensure_ascii=False, default=str)
        self.logger.info(f"批量处理完成！总耗时: {processing_time:.1f} 秒"); self.logger.info(f"最终报告: {report_file}")
        return final_summary


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='批量处理发育阶段接触衰减曲线分析')
    
    # 必需参数
    parser.add_argument('-i', '--input', required=True, help='输入数据根目录 (包含各个stage子目录)')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    # *** 新增的必需参数 ***
    parser.add_argument('-r', '--resolution', required=True, help='指定要处理的分辨率目录，例如 "100K"')
    
    # 可选参数
    parser.add_argument('-s', '--stages', nargs='+', help='指定要处理的发育阶段 (默认处理所有发现的阶段)')
    parser.add_argument('-n', '--max-files', type=int, help='每个阶段最大处理文件数量')
    parser.add_argument('-w', '--max-workers', type=int, default=4, help='并行处理的最大worker数量')
    parser.add_argument('--no-parallel', action='store_true', help='禁用并行处理')
    
    # 分析参数
    parser.add_argument('-c', '--chrom', help='分析指定染色体')
    parser.add_argument('-d', '--max-distance', type=int, help='最大分析距离')
    parser.add_argument('--no-plots', action='store_true', help='不生成图表')
    
    args = parser.parse_args()
    
    # 创建处理器，并传入 resolution 参数
    processor = StageWiseDecayProcessor(
        data_root_dir=args.input,
        output_dir=args.output,
        resolution=args.resolution, # <--- 将 resolution 传递给处理器
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
    print(f"📊 处理了 {results['stages_processed']} 个发育阶段 @ {args.resolution} 分辨率")
    print(f"📁 成功分析 {results['total_successful']}/{results['total_files']} 个文件")
    print(f"⏱️  总耗时: {results['processing_time_seconds']:.1f} 秒")
    print(f"📈 成功率: {results['overall_success_rate']*100:.1f}%")
    print(f"📂 输出目录: {results['output_directory']}")

if __name__ == "__main__":
    main()