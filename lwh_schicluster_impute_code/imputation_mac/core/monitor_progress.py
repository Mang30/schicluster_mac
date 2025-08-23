#!/usr/bin/env python3
"""
M4 MacBook MPS 插补进度监控脚本
实时监控各任务进度、资源使用率和性能统计
"""

import os
import time
import subprocess
import psutil
import json
from datetime import datetime
import argparse

class ProgressMonitor:
    def __init__(self):
        self.output_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage"
        self.log_dir = "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac"
        self.stages = ["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"]
        self.expected_totals = {
            "E70": 557, "E75": 1870, "E80": 559, "E85": 1125,
            "E95": 1108, "EX05": 1128, "EX15": 1122
        }
        
    def get_system_stats(self):
        """获取系统资源使用情况"""
        try:
            cpu_percent = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory()
            
            # 查找Python进程（插补任务）
            python_processes = []
            for proc in psutil.process_iter(['pid', 'name', 'cpu_percent', 'memory_percent']):
                try:
                    if 'python' in proc.info['name'].lower():
                        python_processes.append({
                            'pid': proc.info['pid'],
                            'cpu': proc.info['cpu_percent'],
                            'memory': proc.info['memory_percent']
                        })
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            
            return {
                'cpu_total': cpu_percent,
                'memory_total': memory.percent,
                'memory_used_gb': memory.used / (1024**3),
                'memory_total_gb': memory.total / (1024**3),
                'python_processes': len(python_processes),
                'python_cpu_sum': sum(p['cpu'] for p in python_processes),
                'python_memory_sum': sum(p['memory'] for p in python_processes)
            }
        except Exception as e:
            print(f"获取系统统计失败: {e}")
            return None
    
    def get_progress_stats(self):
        """获取各阶段插补进度"""
        stats = {}
        total_completed = 0
        total_expected = 0
        
        for stage in self.stages:
            stage_dir = os.path.join(self.output_dir, stage)
            completed = 0
            
            if os.path.exists(stage_dir):
                # 统计已完成的矩阵文件
                try:
                    result = subprocess.run(
                        ["find", stage_dir, "-name", "*.npz"],
                        capture_output=True, text=True
                    )
                    if result.stdout.strip():
                        completed = len(result.stdout.strip().split('\n'))
                except:
                    completed = 0
            
            expected = self.expected_totals.get(stage, 0) * 20  # 每细胞20个染色体
            progress = (completed / expected * 100) if expected > 0 else 0
            
            stats[stage] = {
                'completed': completed,
                'expected': expected,
                'progress': progress
            }
            
            total_completed += completed
            total_expected += expected
        
        overall_progress = (total_completed / total_expected * 100) if total_expected > 0 else 0
        
        return {
            'stages': stats,
            'total_completed': total_completed,
            'total_expected': total_expected,
            'overall_progress': overall_progress
        }
    
    def get_task_status(self):
        """检查各任务脚本的运行状态"""
        task_status = {}
        
        for i in range(1, 7):
            log_file = f"task{i}_imputation.log"
            log_path = os.path.join(self.log_dir, log_file)
            
            status = "未开始"
            last_update = "N/A"
            
            if os.path.exists(log_path):
                try:
                    # 读取最后几行日志判断状态
                    with open(log_path, 'r') as f:
                        lines = f.readlines()
                    
                    if lines:
                        last_line = lines[-1].strip()
                        last_update = datetime.fromtimestamp(
                            os.path.getmtime(log_path)
                        ).strftime("%H:%M:%S")
                        
                        if "执行成功" in last_line:
                            status = "已完成"
                        elif "执行失败" in last_line:
                            status = "失败"
                        elif "开始" in last_line or any("完成" in line for line in lines[-5:]):
                            status = "运行中"
                except:
                    status = "无法读取"
            
            task_status[f"Task{i}"] = {
                'status': status,
                'last_update': last_update,
                'log_file': log_file
            }
        
        return task_status
    
    def display_dashboard(self, clear_screen=True):
        """显示实时监控面板"""
        if clear_screen:
            os.system('clear')
        
        print("🚀 M4 MacBook 24GB MPS 插补进度监控")
        print("=" * 80)
        print(f"📅 监控时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print()
        
        # 系统资源
        sys_stats = self.get_system_stats()
        if sys_stats:
            print("💻 系统资源:")
            print(f"   CPU使用率: {sys_stats['cpu_total']:.1f}%")
            print(f"   内存使用: {sys_stats['memory_used_gb']:.1f}GB / {sys_stats['memory_total_gb']:.1f}GB ({sys_stats['memory_total']:.1f}%)")
            print(f"   Python进程: {sys_stats['python_processes']} 个")
            if sys_stats['python_processes'] > 0:
                print(f"   插补任务CPU: {sys_stats['python_cpu_sum']:.1f}%")
                print(f"   插补任务内存: {sys_stats['python_memory_sum']:.1f}%")
        print()
        
        # 任务状态
        task_status = self.get_task_status()
        print("📋 任务状态:")
        for task_name, info in task_status.items():
            status_emoji = {
                "未开始": "⏸️",
                "运行中": "🔄",
                "已完成": "✅",
                "失败": "❌",
                "无法读取": "⚠️"
            }
            emoji = status_emoji.get(info['status'], "❓")
            print(f"   {emoji} {task_name}: {info['status']} (最后更新: {info['last_update']})")
        print()
        
        # 插补进度
        progress_stats = self.get_progress_stats()
        print("📊 插补进度:")
        for stage, stats in progress_stats['stages'].items():
            bar_length = 20
            filled_length = int(bar_length * stats['progress'] / 100)
            bar = "█" * filled_length + "░" * (bar_length - filled_length)
            print(f"   {stage}: [{bar}] {stats['progress']:5.1f}% ({stats['completed']:6d}/{stats['expected']:6d})")
        
        print()
        print(f"📈 总体进度: {progress_stats['overall_progress']:.1f}% ({progress_stats['total_completed']:,}/{progress_stats['total_expected']:,})")
        
        # 预估时间
        if progress_stats['overall_progress'] > 5:  # 至少完成5%才预估
            # 简单的时间预估（基于已完成的比例）
            current_time = time.time()
            # 这里需要记录开始时间，暂时用简化版本
            print(f"⏱️  预估: 根据当前进度预估剩余时间...")
        
        print("=" * 80)
    
    def save_progress_snapshot(self):
        """保存当前进度快照到JSON文件"""
        snapshot = {
            'timestamp': datetime.now().isoformat(),
            'system': self.get_system_stats(),
            'progress': self.get_progress_stats(),
            'tasks': self.get_task_status()
        }
        
        snapshot_file = os.path.join(self.log_dir, "progress_snapshot.json")
        try:
            with open(snapshot_file, 'w') as f:
                json.dump(snapshot, f, indent=2)
        except Exception as e:
            print(f"保存快照失败: {e}")
    
    def monitor_loop(self, interval=30, save_snapshots=False):
        """持续监控循环"""
        try:
            while True:
                self.display_dashboard()
                
                if save_snapshots:
                    self.save_progress_snapshot()
                
                # 检查是否所有任务完成
                task_status = self.get_task_status()
                running_tasks = [t for t, info in task_status.items() if info['status'] == '运行中']
                
                if not running_tasks:
                    progress_stats = self.get_progress_stats()
                    if progress_stats['overall_progress'] > 95:
                        print("\n🎉 所有任务可能已完成!")
                        break
                
                print(f"\n⏱️  下次更新: {interval}秒后 (按 Ctrl+C 退出)")
                time.sleep(interval)
                
        except KeyboardInterrupt:
            print("\n👋 监控已停止")
    
    def generate_summary_report(self):
        """生成汇总报告"""
        print("📄 生成最终汇总报告...")
        
        progress_stats = self.get_progress_stats()
        task_status = self.get_task_status()
        sys_stats = self.get_system_stats()
        
        report_file = os.path.join(self.log_dir, f"imputation_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
        
        with open(report_file, 'w') as f:
            f.write("M4 MacBook 24GB MPS 插补汇总报告\n")
            f.write("=" * 60 + "\n")
            f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("📊 最终进度统计:\n")
            for stage, stats in progress_stats['stages'].items():
                f.write(f"   {stage}: {stats['completed']:6d}/{stats['expected']:6d} ({stats['progress']:5.1f}%)\n")
            
            f.write(f"\n总体完成度: {progress_stats['overall_progress']:.1f}%\n")
            f.write(f"总矩阵文件: {progress_stats['total_completed']:,} 个\n\n")
            
            f.write("📋 任务执行状态:\n")
            for task_name, info in task_status.items():
                f.write(f"   {task_name}: {info['status']}\n")
            
            if sys_stats:
                f.write(f"\n💻 系统资源使用:\n")
                f.write(f"   内存使用: {sys_stats['memory_used_gb']:.1f}GB / {sys_stats['memory_total_gb']:.1f}GB\n")
                f.write(f"   Python进程: {sys_stats['python_processes']} 个\n")
        
        print(f"✅ 汇总报告已保存: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='M4 MacBook MPS 插补进度监控')
    parser.add_argument('--interval', type=int, default=30, help='监控间隔(秒)')
    parser.add_argument('--once', action='store_true', help='只显示一次状态')
    parser.add_argument('--save-snapshots', action='store_true', help='保存进度快照')
    parser.add_argument('--summary', action='store_true', help='生成汇总报告')
    
    args = parser.parse_args()
    
    monitor = ProgressMonitor()
    
    if args.summary:
        monitor.generate_summary_report()
    elif args.once:
        monitor.display_dashboard()
    else:
        monitor.monitor_loop(args.interval, args.save_snapshots)

if __name__ == "__main__":
    main()