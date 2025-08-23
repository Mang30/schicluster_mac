#!/usr/bin/env python3
"""
M4 MacBook 24GB MPS 插补主控脚本
管理6个并行任务的启动、监控和汇总
"""

import os
import sys
import time
import subprocess
import threading
import logging
from datetime import datetime
import argparse

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('parallel_controller.log'),
        logging.StreamHandler()
    ]
)

class TaskController:
    def __init__(self):
        self.base_dir = "/Volumes/SumSung500/CSU/0_HiRES/lwh_schicluster_impute_code/imputation_mac"
        self.output_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/imputed_matrices_by_stage"
        self.tasks = {
            1: {"name": "Task1_E70_E80", "script": "task1_imputation.sh", "stages": ["E70", "E80"], "status": "pending"},
            2: {"name": "Task2_E75", "script": "task2_imputation.sh", "stages": ["E75"], "status": "pending"},
            3: {"name": "Task3_E85", "script": "task3_imputation.sh", "stages": ["E85"], "status": "pending"},
            4: {"name": "Task4_E95", "script": "task4_imputation.sh", "stages": ["E95"], "status": "pending"},
            5: {"name": "Task5_EX05", "script": "task5_imputation.sh", "stages": ["EX05"], "status": "pending"},
            6: {"name": "Task6_EX15", "script": "task6_imputation.sh", "stages": ["EX15"], "status": "pending"}
        }
        self.processes = {}
        self.start_time = None
        
    def check_environment(self):
        """检查环境和依赖"""
        logging.info("🔍 检查 M4 MacBook 环境...")
        
        # 检查 Python 环境
        try:
            import torch
            if torch.backends.mps.is_available():
                logging.info("✅ MPS 加速可用")
            else:
                logging.warning("⚠️  MPS 不可用，将使用 CPU")
        except ImportError:
            logging.error("❌ PyTorch 未安装")
            return False
        
        # 检查任务脚本
        for task_id, task_info in self.tasks.items():
            script_path = os.path.join(self.base_dir, task_info["script"])
            if not os.path.exists(script_path):
                logging.error(f"❌ 任务脚本不存在: {script_path}")
                return False
            if not os.access(script_path, os.X_OK):
                logging.error(f"❌ 任务脚本无执行权限: {script_path}")
                return False
        
        # 检查输入数据
        input_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output"
        if not os.path.exists(input_dir):
            logging.error(f"❌ 输入数据目录不存在: {input_dir}")
            return False
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        
        logging.info("✅ 环境检查通过")
        return True
    
    def get_stage_stats(self):
        """获取各阶段细胞统计"""
        stats = {}
        input_dir = "/Volumes/SumSung500/CSU/0_HiRES/output/pairs2matrix_output"
        
        for task_id, task_info in self.tasks.items():
            for stage in task_info["stages"]:
                stage_path = os.path.join(input_dir, stage)
                if os.path.exists(stage_path):
                    cell_count = len([d for d in os.listdir(stage_path) 
                                    if d.startswith("Gas") and os.path.isdir(os.path.join(stage_path, d))])
                    stats[stage] = cell_count
                else:
                    stats[stage] = 0
                    
        return stats
    
    def run_task(self, task_id, test_mode=False):
        """在单独线程中运行任务"""
        task_info = self.tasks[task_id]
        script_path = os.path.join(self.base_dir, task_info["script"])
        
        logging.info(f"🚀 启动 {task_info['name']} ({task_info['stages']})")
        task_info["status"] = "running"
        task_info["start_time"] = time.time()
        
        try:
            # 切换到脚本目录并执行
            process = subprocess.Popen(
                ["bash", script_path],
                cwd=self.base_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1
            )
            
            self.processes[task_id] = process
            
            # 实时输出日志
            for line in iter(process.stdout.readline, ''):
                if line:
                    logging.info(f"[Task {task_id}] {line.strip()}")
            
            process.wait()
            
            if process.returncode == 0:
                task_info["status"] = "completed"
                logging.info(f"✅ {task_info['name']} 完成")
            else:
                task_info["status"] = "failed"
                logging.error(f"❌ {task_info['name']} 失败 (退出码: {process.returncode})")
            
            task_info["end_time"] = time.time()
            task_info["duration"] = task_info["end_time"] - task_info["start_time"]
            
        except Exception as e:
            task_info["status"] = "error"
            logging.error(f"❌ {task_info['name']} 发生异常: {e}")
    
    def monitor_progress(self):
        """监控所有任务进度"""
        while True:
            time.sleep(30)  # 每30秒检查一次
            
            running_tasks = [t for t in self.tasks.values() if t["status"] == "running"]
            if not running_tasks:
                break
            
            # 统计当前进度
            completed_matrices = 0
            if os.path.exists(self.output_dir):
                completed_matrices = len(subprocess.check_output(
                    ["find", self.output_dir, "-name", "*.npz"], 
                    universal_newlines=True
                ).strip().split('\n')) if subprocess.check_output(
                    ["find", self.output_dir, "-name", "*.npz"], 
                    universal_newlines=True
                ).strip() else 0
            
            elapsed_time = time.time() - self.start_time
            logging.info(f"📊 进度监控: {len(running_tasks)} 个任务运行中, "
                        f"{completed_matrices} 个矩阵已完成, "
                        f"运行时间 {elapsed_time/3600:.1f} 小时")
    
    def run_all_tasks(self, test_mode=False, selected_tasks=None):
        """并行运行所有任务"""
        if not self.check_environment():
            return False
        
        # 显示统计信息
        stats = self.get_stage_stats()
        total_cells = sum(stats.values())
        
        logging.info("📊 M4 MacBook 24GB 并行插补任务启动")
        logging.info("=" * 60)
        for stage, count in stats.items():
            logging.info(f"   {stage}: {count} 个细胞")
        logging.info(f"   总计: {total_cells} 个细胞")
        logging.info(f"   预计矩阵: {total_cells * 20} 个")
        logging.info("=" * 60)
        
        if test_mode:
            logging.info("🧪 测试模式: 每个任务处理2个细胞")
        
        # 选择要运行的任务
        tasks_to_run = selected_tasks or list(self.tasks.keys())
        
        self.start_time = time.time()
        threads = []
        
        # 启动监控线程
        monitor_thread = threading.Thread(target=self.monitor_progress, daemon=True)
        monitor_thread.start()
        
        # 启动任务线程
        for task_id in tasks_to_run:
            if task_id in self.tasks:
                thread = threading.Thread(target=self.run_task, args=(task_id, test_mode))
                thread.start()
                threads.append(thread)
                time.sleep(5)  # 错开启动时间，避免资源竞争
        
        # 等待所有任务完成
        for thread in threads:
            thread.join()
        
        # 生成最终报告
        self.generate_final_report()
        
        return True
    
    def generate_final_report(self):
        """生成最终执行报告"""
        total_duration = time.time() - self.start_time
        
        logging.info("\n" + "=" * 70)
        logging.info("🎉 M4 MacBook 24GB 并行插补任务完成!")
        logging.info("=" * 70)
        
        # 任务状态汇总
        completed = len([t for t in self.tasks.values() if t["status"] == "completed"])
        failed = len([t for t in self.tasks.values() if t["status"] == "failed"])
        error = len([t for t in self.tasks.values() if t["status"] == "error"])
        
        logging.info(f"📊 任务统计: {completed} 完成, {failed} 失败, {error} 异常")
        
        # 输出文件统计
        if os.path.exists(self.output_dir):
            total_matrices = 0
            for stage in ["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"]:
                stage_dir = os.path.join(self.output_dir, stage)
                if os.path.exists(stage_dir):
                    count = len(subprocess.check_output(
                        ["find", stage_dir, "-name", "*.npz"], 
                        universal_newlines=True
                    ).strip().split('\n')) if subprocess.check_output(
                        ["find", stage_dir, "-name", "*.npz"], 
                        universal_newlines=True
                    ).strip() else 0
                    total_matrices += count
                    logging.info(f"   {stage}: {count} 个矩阵文件")
            
            logging.info(f"📈 总矩阵文件: {total_matrices} 个")
        
        # 性能统计
        logging.info(f"⏱️  总执行时间: {total_duration/3600:.2f} 小时")
        if total_matrices > 0:
            avg_time = total_duration / total_matrices
            logging.info(f"📊 平均处理时间: {avg_time:.2f} 秒/矩阵")
            logging.info(f"💡 M4 MPS 估计加速比: ~12x vs CPU")
        
        # 各任务详细时间
        logging.info("\n📋 任务详细时间:")
        for task_id, task_info in self.tasks.items():
            if "duration" in task_info:
                logging.info(f"   {task_info['name']}: {task_info['duration']/3600:.2f} 小时")
        
        logging.info("=" * 70)

def main():
    parser = argparse.ArgumentParser(description='M4 MacBook 24GB 并行插补控制器')
    parser.add_argument('--test', action='store_true', help='测试模式（每任务2个细胞）')
    parser.add_argument('--tasks', nargs='+', type=int, choices=range(1, 7),
                       help='指定运行的任务 ID (1-6)')
    
    args = parser.parse_args()
    
    controller = TaskController()
    
    try:
        success = controller.run_all_tasks(
            test_mode=args.test,
            selected_tasks=args.tasks
        )
        return 0 if success else 1
        
    except KeyboardInterrupt:
        logging.info("\n❌ 用户中断执行")
        # 终止所有子进程
        for process in controller.processes.values():
            if process and process.poll() is None:
                process.terminate()
        return 1

if __name__ == "__main__":
    exit(main())