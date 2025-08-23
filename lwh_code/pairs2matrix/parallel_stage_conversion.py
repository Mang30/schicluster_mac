#!/usr/bin/env python3
"""
Parallel stage-specific conversion script for pairs files
"""

import os
import glob
import re
import argparse
from pathlib import Path
import subprocess

def extract_stage_and_cell_id(filepath):
    """Extract stage and cell ID from filepath"""
    filename = os.path.basename(filepath)
    
    # Gas* pattern: Gas[a-z]E[stage][cell]
    gas_match = re.search(r'Gas[a-z]E(\d{2})(\d{4})', filename)
    if gas_match:
        stage = f"E{gas_match.group(1)}"
        cell_number = gas_match.group(2)
        prefix_match = re.search(r'(Gas[a-z]E)', filename)
        if prefix_match:
            cell_id = f"{prefix_match.group(1)}{gas_match.group(1)}{cell_number}"
            return stage, cell_id
    
    # Org* pattern: Org[a-z]E[stage][cell] 
    org_match = re.search(r'Org[a-z]E(\d{2})(\d{4})', filename)
    if org_match:
        stage = f"E{org_match.group(1)}"
        cell_number = org_match.group(2)
        prefix_match = re.search(r'(Org[a-z]E)', filename)
        if prefix_match:
            cell_id = f"{prefix_match.group(1)}{org_match.group(1)}{cell_number}"
            return stage, cell_id
    
    # EX pattern: Org[a-z]EX[stage][cell]
    ex_match = re.search(r'Org[a-z]EX(\d{2})(\d{4})', filename)
    if ex_match:
        stage = f"EX{ex_match.group(1)}"
        cell_number = ex_match.group(2)
        prefix_match = re.search(r'(Org[a-z]EX)', filename)
        if prefix_match:
            cell_id = f"{prefix_match.group(1)}{ex_match.group(1)}{cell_number}"
            return stage, cell_id
    
    return None, None

def convert_stage_data(input_dir, output_base_dir, chrom_sizes, target_stages, resolution=100000):
    """
    Convert pairs files for specific stages
    """
    
    # Find all pairs files
    pairs_files = glob.glob(os.path.join(input_dir, "*.pairs.gz"))
    
    # Group files by stage
    stage_files = {}
    for pairs_file in pairs_files:
        stage, cell_id = extract_stage_and_cell_id(pairs_file)
        if stage and cell_id and stage in target_stages:
            if stage not in stage_files:
                stage_files[stage] = []
            stage_files[stage].append((pairs_file, cell_id))
    
    print(f"Processing stages: {target_stages}")
    for stage in target_stages:
        if stage in stage_files:
            print(f"  {stage}: {len(stage_files[stage])} cells")
        else:
            print(f"  {stage}: 0 cells found")
    
    # Process each target stage
    total_processed = 0
    for stage in target_stages:
        if stage not in stage_files:
            continue
            
        files_for_stage = stage_files[stage]
        print(f"\\nProcessing {len(files_for_stage)} files for {stage}")
        
        # Create stage-specific output directory
        stage_output_dir = os.path.join(output_base_dir, stage)
        os.makedirs(stage_output_dir, exist_ok=True)
        
        # Convert each file in this stage
        for i, (pairs_file, cell_id) in enumerate(files_for_stage):
            cell_output_dir = os.path.join(stage_output_dir, cell_id)
            
            print(f"  [{i+1}/{len(files_for_stage)}] Converting {cell_id}...")
            
            # Run conversion script
            cmd = [
                'python', 'lwh_code/pairs2matrix/convert_pairs_to_matrix.py',
                '--pairs_file', pairs_file,
                '--output_dir', cell_output_dir,
                '--cell_id', cell_id,
                '--resolution', str(resolution),
                '--chrom_sizes', chrom_sizes
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    print(f"    ✓ Successfully converted")
                    total_processed += 1
                else:
                    print(f"    ✗ Failed to convert")
                    print(f"    Error: {result.stderr}")
            except Exception as e:
                print(f"    ✗ Exception: {e}")
        
        print(f"  Stage {stage} completed")
    
    print(f"\\n=== Total processed: {total_processed} cells ===")

def main():
    parser = argparse.ArgumentParser(description='Parallel stage-specific conversion')
    parser.add_argument('--input_dir', required=True, help='Directory containing .pairs.gz files')
    parser.add_argument('--output_dir', default='output/pairs2matrix_output', help='Output base directory')
    parser.add_argument('--chrom_sizes', required=True, help='Chromosome sizes file')
    parser.add_argument('--stages', nargs='+', required=True, help='Target stages to process')
    parser.add_argument('--resolution', type=int, default=100000, help='Resolution in bp')
    
    args = parser.parse_args()
    
    convert_stage_data(
        args.input_dir,
        args.output_dir, 
        args.chrom_sizes,
        args.stages,
        args.resolution
    )

if __name__ == "__main__":
    main()