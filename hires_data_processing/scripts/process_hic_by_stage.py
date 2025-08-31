#!/usr/bin/env python3
"""
Process Hi-C data using schicluster for each developmental stage
"""

import os
import subprocess
import argparse
from pathlib import Path
import logging
from datetime import datetime

def setup_logging(log_dir):
    """Setup logging configuration"""
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"schicluster_processing_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def run_schicluster_prepare(stage, contact_table, output_dir, chrom_size_path, 
                           resolution=100000, batch_size=1536, cpu_per_job=30,
                           pad=1, chr1=1, pos1=2, chr2=5, pos2=6, 
                           output_dist=500000000, window_size=500000000, 
                           step_size=500000000, logger=None):
    """Run schicluster prepare-impute for a specific stage"""
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Prepare command
    cmd = [
        "hicluster", "prepare-impute",
        "--cell_table", contact_table,
        "--batch_size", str(batch_size),
        "--pad", str(pad),
        "--cpu_per_job", str(cpu_per_job),
        "--chr1", str(chr1),
        "--pos1", str(pos1),
        "--chr2", str(chr2),
        "--pos2", str(pos2),
        "--output_dir", output_dir,
        "--chrom_size_path", chrom_size_path,
        "--output_dist", str(output_dist),
        "--window_size", str(window_size),
        "--step_size", str(step_size),
        "--resolution", str(resolution)
    ]
    
    logger.info(f"Step 1: Running prepare-impute for stage {stage}")
    logger.info(f"Command: {' '.join(cmd)}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Contact table: {contact_table}")
    
    try:
        # Run the command (compatible with older Python versions)
        result = subprocess.run(
            cmd, 
            cwd=os.path.dirname(output_dir),
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            universal_newlines=True, 
            timeout=3600  # 1 hour timeout for prepare step
        )
        
        if result.returncode == 0:
            logger.info(f"Successfully completed prepare-impute for stage {stage}")
            logger.info(f"STDOUT: {result.stdout}")
            return True
        else:
            logger.error(f"Error in prepare-impute for stage {stage}")
            logger.error(f"Return code: {result.returncode}")
            logger.error(f"STDERR: {result.stderr}")
            logger.error(f"STDOUT: {result.stdout}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error(f"Timeout expired in prepare-impute for stage {stage}")
        return False
    except Exception as e:
        logger.error(f"Exception in prepare-impute for stage {stage}: {str(e)}")
        return False

def run_snakemake_imputation(stage, output_dir, cpu_per_job=30, logger=None):
    """Run the snakemake workflow for imputation"""
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    # Find the snakemake command file
    snakemake_cmd_file = os.path.join(output_dir, "snakemake_cmd.txt")
    if not os.path.exists(snakemake_cmd_file):
        logger.error(f"Snakemake command file not found: {snakemake_cmd_file}")
        return False
    
    # Read the snakemake command
    with open(snakemake_cmd_file, 'r') as f:
        snakemake_cmd = f.read().strip()
    
    logger.info(f"Step 2: Running snakemake imputation for stage {stage}")
    logger.info(f"Command: {snakemake_cmd}")
    
    try:
        # Parse and run the snakemake command (compatible with older Python versions)
        cmd_parts = snakemake_cmd.split()
        result = subprocess.run(
            cmd_parts,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=86400  # 24 hour timeout for imputation
        )
        
        if result.returncode == 0:
            logger.info(f"Successfully completed snakemake imputation for stage {stage}")
            logger.info(f"STDOUT: {result.stdout}")
            return True
        else:
            logger.error(f"Error in snakemake imputation for stage {stage}")
            logger.error(f"Return code: {result.returncode}")
            logger.error(f"STDERR: {result.stderr}")
            logger.error(f"STDOUT: {result.stdout}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error(f"Timeout expired in snakemake imputation for stage {stage}")
        return False
    except Exception as e:
        logger.error(f"Exception in snakemake imputation for stage {stage}: {str(e)}")
        return False

def run_schicluster_full_pipeline(stage, contact_table, output_dir, chrom_size_path, 
                                 resolution=100000, batch_size=1536, cpu_per_job=30,
                                 pad=1, chr1=1, pos1=2, chr2=3, pos2=4, 
                                 output_dist=500000000, window_size=500000000, 
                                 step_size=500000000, logger=None):
    """Run the complete schicluster pipeline: prepare + imputation"""
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f"Starting complete schicluster pipeline for stage {stage}")
    
    # Step 1: Prepare imputation
    success = run_schicluster_prepare(
        stage, contact_table, output_dir, chrom_size_path,
        resolution, batch_size, cpu_per_job, pad, chr1, pos1, chr2, pos2,
        output_dist, window_size, step_size, logger
    )
    
    if not success:
        logger.error(f"Prepare step failed for stage {stage}")
        return False
    
    # Step 2: Run snakemake imputation
    success = run_snakemake_imputation(stage, output_dir, cpu_per_job, logger)
    
    if not success:
        logger.error(f"Imputation step failed for stage {stage}")
        return False
    
    logger.info(f"Complete pipeline successful for stage {stage}")
    return True

def get_cell_count(contact_table):
    """Get number of cells in contact table"""
    try:
        with open(contact_table, 'r') as f:
            return sum(1 for line in f if line.strip())
    except:
        return 0

def main():
    parser = argparse.ArgumentParser(description="Process Hi-C data using schicluster by developmental stage")
    parser.add_argument("--stages", nargs="+", default=["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"],
                       help="Developmental stages to process")
    parser.add_argument("--resolution", type=int, default=100000,
                       help="Resolution for analysis (default: 100000)")
    parser.add_argument("--batch_size", type=int, default=5000,
                       help="Batch size (default: 5000)")
    parser.add_argument("--cpu_per_job", type=int, default=7,
                       help="CPUs per job (default: 30)")
    parser.add_argument("--specific_stage", type=str, default=None,
                       help="Process only a specific stage")
    parser.add_argument("--chrom_size_file", type=str, default=None,
                       help="Chromosome size file to use (default: mm10_chrom_sizes_with_chrY.txt)")
    
    args = parser.parse_args()
    
    # Use relative paths from the script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    base_dir = project_dir
    contact_table_dir = os.path.join(base_dir, "contact_tables")
    output_base_dir = os.path.join(base_dir, "outputs")
    
    # Determine chromosome size file to use
    if args.chrom_size_file:
        chrom_size_path = os.path.join(base_dir, args.chrom_size_file)
    else:
        chrom_size_path = os.path.join(base_dir, "mm10_chrom_sizes_with_chrY.txt")
    
    log_dir = os.path.join(base_dir, "logs")
    
    # Setup logging
    logger = setup_logging(log_dir)
    
    # Verify chromosome size file exists
    if not os.path.exists(chrom_size_path):
        logger.error(f"Chromosome size file not found: {chrom_size_path}")
        return
    
    # Determine which stages to process
    stages_to_process = [args.specific_stage] if args.specific_stage else args.stages
    
    logger.info(f"Starting schicluster processing for stages: {stages_to_process}")
    logger.info(f"Parameters: resolution={args.resolution}, batch_size={args.batch_size}, cpu_per_job={args.cpu_per_job}")
    
    results = {}
    
    for stage in stages_to_process:
        logger.info(f"\n{'='*50}")
        logger.info(f"Processing stage: {stage}")
        logger.info(f"{'='*50}")
        
        # Paths for this stage
        # Handle special case for length-specific contact tables
        if stage.endswith("length"):
            contact_table = os.path.join(contact_table_dir, f"{stage}_chromosome_contact_table.tsv")
        else:
            contact_table = os.path.join(contact_table_dir, f"contact_table_{stage}.tsv")
            
        output_dir = os.path.join(output_base_dir, stage, "impute", f"{args.resolution//1000}K")
        
        # Check if contact table exists
        if not os.path.exists(contact_table):
            logger.error(f"Contact table not found: {contact_table}")
            results[stage] = "Contact table missing"
            continue
            
        # Get cell count
        cell_count = get_cell_count(contact_table)
        logger.info(f"Processing {cell_count} cells for stage {stage}")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Run complete schicluster pipeline
        success = run_schicluster_full_pipeline(
            stage=stage,
            contact_table=contact_table,
            output_dir=output_dir,
            chrom_size_path=chrom_size_path,
            resolution=args.resolution,
            batch_size=args.batch_size,
            cpu_per_job=args.cpu_per_job,
            logger=logger
        )
        
        results[stage] = "Success" if success else "Failed"
        
        logger.info(f"Stage {stage} processing result: {results[stage]}")
    
    # Summary
    logger.info(f"\n{'='*50}")
    logger.info("PROCESSING SUMMARY")
    logger.info(f"{'='*50}")
    
    for stage, result in results.items():
        logger.info(f"{stage}: {result}")
    
    successful_stages = sum(1 for result in results.values() if result == "Success")
    logger.info(f"\nSuccessfully processed {successful_stages}/{len(results)} stages")

if __name__ == "__main__":
    main()