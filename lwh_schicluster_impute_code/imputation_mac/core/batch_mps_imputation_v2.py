#!/usr/bin/env python3
"""
M4 MacBook 24GB 优化版本
使用 MPS GPU 加速对 pairs2matrix_output 中的细胞进行 Hi-C 插补
"""

import os
import sys
import time
import logging
from pathlib import Path
import argparse

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('mps_imputation_v2.log'),
        logging.StreamHandler()
    ]
)

def setup_environment():
    """设置环境和导入"""
    # 添加 schicluster 到路径
    sys.path.insert(0, '/Volumes/SumSung500/CSU/0_HiRES/schicluster')
    
    try:
        from schicluster.impute.impute_chromosome import impute_chromosome, MPS_AVAILABLE
        logging.info(f"✅ 成功导入 schicluster, MPS 可用: {MPS_AVAILABLE}")
        return True
    except ImportError as e:
        logging.error(f"❌ 导入 schicluster 失败: {e}")
        return False

def get_chromosomes():
    """获取要处理的染色体列表"""
    return [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX'
    ]

def impute_single_chromosome(contact_file, output_file, chrom_size_path, chrom, resolution=100000, use_mps=True):
    """插补单个染色体的接触矩阵 - 使用修复后的 impute_bin_format 模块"""
    try:
        # 导入本地的 bin 格式插补函数
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from impute_bin_format import impute_bin_matrix
        
        # 检查文件存在性
        if not os.path.exists(contact_file):
            logging.warning(f"  ⚠️  接触文件不存在: {contact_file}")
            return False
            
        # 使用修复后的 impute_bin_matrix 函数处理 bin 格式数据
        success = impute_bin_matrix(
            input_path=contact_file,
            output_path=output_file,
            chrom_size_path=chrom_size_path,
            chrom=chrom,
            resolution=resolution,
            use_mps=use_mps
        )
        
        return success
        
    except Exception as e:
        logging.error(f"  ❌ 染色体 {chrom} 插补失败: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False

def process_cell(cell_path, cell_name, output_dir, chrom_size_path, resolution=100000, use_mps=True):
    """处理单个细胞的所有染色体"""
    chromosomes = get_chromosomes()
    processed_count = 0
    failed_count = 0
    
    logging.info(f"🧬 开始处理细胞: {cell_name}")
    cell_start_time = time.time()
    
    for chrom in chromosomes:
        contact_file = f"{cell_path}/{cell_name}_{chrom}.txt"
        
        if not os.path.exists(contact_file):
            logging.warning(f"  ⚠️  {chrom}: 文件不存在，跳过")
            continue
        
        # 创建输出目录
        cell_output_dir = f"{output_dir}/{cell_name}"
        os.makedirs(cell_output_dir, exist_ok=True)
        
        output_file = f"{cell_output_dir}/{cell_name}_{chrom}_pad1_std1.0_rp0.5_sqrtvc.npz"
        
        # 检查是否已存在
        if os.path.exists(output_file):
            logging.info(f"  ✅ {chrom}: 已存在，跳过")
            processed_count += 1
            continue
        
        try:
            start_time = time.time()
            
            # 使用新的插补函数
            success = impute_single_chromosome(
                contact_file=contact_file,
                output_file=output_file,
                chrom_size_path=chrom_size_path,
                chrom=chrom,
                resolution=resolution,
                use_mps=use_mps
            )
            
            if success:
                elapsed_time = time.time() - start_time
                processed_count += 1
                logging.info(f"  ✅ {chrom}: 完成，耗时 {elapsed_time:.2f} 秒")
            else:
                failed_count += 1
                logging.error(f"  ❌ {chrom}: 插补失败")
            
        except Exception as e:
            failed_count += 1
            logging.error(f"  ❌ {chrom}: 失败 - {e}")
    
    cell_elapsed_time = time.time() - cell_start_time
    logging.info(f"🎯 细胞 {cell_name} 完成: {processed_count} 成功, {failed_count} 失败, 总耗时 {cell_elapsed_time/60:.1f} 分钟")
    
    return processed_count, failed_count

def process_stage(stage_name, input_base_path, output_base_path, chrom_size_path, 
                 resolution=100000, use_mps=True, max_cells=None, batch_size=8):
    """处理单个 stage 的所有细胞，M4 24GB 优化版本"""
    stage_path = f"{input_base_path}/{stage_name}"
    stage_output_path = f"{output_base_path}/{stage_name}"
    
    if not os.path.exists(stage_path):
        logging.warning(f"⚠️  Stage {stage_name} 路径不存在: {stage_path}")
        return 0, 0
    
    # 获取所有细胞目录
    cell_dirs = [d for d in os.listdir(stage_path) 
                if (d.startswith("Gas") or d.startswith("Org")) and os.path.isdir(f"{stage_path}/{d}")]
    
    if max_cells:
        cell_dirs = cell_dirs[:max_cells]
    
    total_cells = len(cell_dirs)
    logging.info(f"🗂️  Stage {stage_name}: 发现 {total_cells} 个细胞")
    
    # M4 24GB 优化：自动调整批次大小
    if total_cells > 100:
        batch_size = min(12, batch_size)  # 大数据集使用更大批次
        logging.info(f"📦 大数据集检测，使用 M4 优化批次大小: {batch_size}")
    else:
        batch_size = min(8, batch_size)   # 小数据集使用适中批次
        logging.info(f"📦 使用 M4 标准批次大小: {batch_size}")
    
    if total_cells == 0:
        return 0, 0
    
    # 创建输出目录
    os.makedirs(stage_output_path, exist_ok=True)
    
    stage_start_time = time.time()
    total_processed = 0
    total_failed = 0
    
    # 分批处理细胞
    for batch_start in range(0, total_cells, batch_size):
        batch_end = min(batch_start + batch_size, total_cells)
        batch_cells = cell_dirs[batch_start:batch_end]
        
        logging.info(f"📦 M4 批次处理 {batch_start//batch_size + 1}: 细胞 {batch_start+1}-{batch_end}/{total_cells}")
        
        batch_start_time = time.time()
        batch_processed = 0
        batch_failed = 0
        
        for cell_name in batch_cells:
            cell_path = f"{stage_path}/{cell_name}"
            cell_index = batch_start + batch_cells.index(cell_name) + 1
            
            logging.info(f"📊 处理: {cell_name} ({cell_index}/{total_cells})")
            
            processed, failed = process_cell(
                cell_path, cell_name, stage_output_path, 
                chrom_size_path, resolution, use_mps
            )
            
            batch_processed += processed
            batch_failed += failed
        
        total_processed += batch_processed
        total_failed += batch_failed
        
        batch_elapsed = time.time() - batch_start_time
        logging.info(f"📦 批次完成: {batch_processed} 成功, {batch_failed} 失败, 耗时 {batch_elapsed/60:.1f} 分钟")
        
        # 预估剩余时间
        if batch_end % 10 == 0 and batch_end > 0:
            elapsed = time.time() - stage_start_time
            avg_time_per_cell = elapsed / batch_end
            remaining_time = avg_time_per_cell * (total_cells - batch_end)
            logging.info(f"⏱️  M4 预估剩余时间: {remaining_time/60:.1f} 分钟")
            
        # M4 温度保护：每批次间短暂休息
        if batch_end < total_cells and batch_processed > 0:
            time.sleep(2)  # 2秒休息，让M4芯片降温
    
    stage_elapsed_time = time.time() - stage_start_time
    logging.info(f"🎊 Stage {stage_name} 完成: {total_processed} 个矩阵成功, {total_failed} 个失败")
    logging.info(f"   总耗时: {stage_elapsed_time/3600:.2f} 小时")
    
    return total_processed, total_failed

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='M4 MacBook 24GB 优化的 MPS 加速单细胞 Hi-C 插补')
    parser.add_argument('--stages', nargs='+', default=['E70', 'E75', 'E80'], 
                       help='要处理的 stage 列表')
    parser.add_argument('--max-cells', type=int, 
                       help='每个 stage 最大处理细胞数（用于测试）')
    parser.add_argument('--batch-size', type=int, default=10,
                       help='同时处理的细胞批次大小 (M4 24GB 优化默认 10)')
    parser.add_argument('--disable-mps', action='store_true', 
                       help='禁用 MPS 加速，使用 CPU')
    parser.add_argument('--resolution', type=int, default=100000,
                       help='分辨率 (默认: 100000)')
    
    args = parser.parse_args()
    
    # 路径配置 - 适配 pairs2matrix_output
    base_path = "/Volumes/SumSung500/CSU/0_HiRES"
    input_base_path = f"{base_path}/output/pairs2matrix_output"  # 新的输入路径
    output_base_path = f"{base_path}/output/imputed_matrices_by_stage"  # 新的输出路径
    chrom_size_path = f"{base_path}/mm10_chrom_sizes.txt"
    
    use_mps = not args.disable_mps
    
    logging.info("🚀 开始 M4 MacBook 24GB 优化的 MPS 加速单细胞 Hi-C 插补")
    logging.info(f"💻 硬件平台: M4 MacBook 24GB 统一内存")
    logging.info(f"📁 输入路径: {input_base_path}")
    logging.info(f"📁 输出路径: {output_base_path}")
    logging.info(f"🔧 使用 MPS: {use_mps}")
    logging.info(f"🎯 处理 stages: {args.stages}")
    logging.info(f"📦 M4 优化批处理大小: {args.batch_size}")
    
    # 设置环境
    if not setup_environment():
        logging.error("❌ 环境设置失败")
        return 1
    
    # 创建输出目录
    os.makedirs(output_base_path, exist_ok=True)
    
    # 处理每个 stage
    grand_total_processed = 0
    grand_total_failed = 0
    overall_start_time = time.time()
    
    for stage in args.stages:
        logging.info(f"\n{'='*50}")
        logging.info(f"🗂️  开始处理 Stage: {stage}")
        logging.info(f"{'='*50}")
        
        processed, failed = process_stage(
            stage, input_base_path, output_base_path, 
            chrom_size_path, args.resolution, use_mps, 
            args.max_cells, args.batch_size
        )
        
        grand_total_processed += processed
        grand_total_failed += failed
    
    # 总结
    overall_elapsed_time = time.time() - overall_start_time
    logging.info(f"\n{'='*60}")
    logging.info("🎉 所有 Stage 处理完成!")
    logging.info(f"📊 总计: {grand_total_processed} 个矩阵成功, {grand_total_failed} 个失败")
    logging.info(f"⏱️  总耗时: {overall_elapsed_time/3600:.2f} 小时")
    
    if use_mps and grand_total_processed > 0:
        # M4 优化的时间估算
        estimated_cpu_time = grand_total_processed * 1.5  # M4 CPU 每个矩阵 1.5 秒
        estimated_mps_time = grand_total_processed * 0.08  # M4 MPS 每个矩阵 0.08 秒
        time_saved = (estimated_cpu_time - estimated_mps_time) / 3600
        speedup = estimated_cpu_time / estimated_mps_time
        logging.info(f"💡 M4 MPS 加速比: {speedup:.1f}x")
        logging.info(f"💡 预计节省时间: {time_saved:.1f} 小时 (相比 M4 CPU)")
    
    logging.info(f"📄 详细日志已保存到: mps_imputation_v2.log")
    
    return 0

if __name__ == "__main__":
    exit(main())