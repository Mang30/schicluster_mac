#!/usr/bin/env python3
"""
Stage-aware conversion script for pairs files to schicluster format
Organizes data by developmental stages
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
    
    # Pattern examples:
    # GSM6998595_GasaE751001.pairs.gz -> E75
    # GSM7001674_OrgaE952002.pairs.gz -> E95  
    # GSM7002031_OrgaEX151001.pairs.gz -> EX15
    
    # Try Gas* pattern first: Gas[a-z]E[stage][cell]
    gas_match = re.search(r'Gas[a-z]E(\d{2})(\d{4})', filename)
    if gas_match:
        stage = f"E{gas_match.group(1)}"
        cell_number = gas_match.group(2)
        # Extract full prefix
        prefix_match = re.search(r'(Gas[a-z]E)', filename)
        if prefix_match:
            cell_id = f"{prefix_match.group(1)}{gas_match.group(1)}{cell_number}"
            return stage, cell_id
    
    # Try Org* pattern: Org[a-z]E[stage][cell] 
    org_match = re.search(r'Org[a-z]E(\d{2})(\d{4})', filename)
    if org_match:
        stage = f"E{org_match.group(1)}"
        cell_number = org_match.group(2)
        prefix_match = re.search(r'(Org[a-z]E)', filename)
        if prefix_match:
            cell_id = f"{prefix_match.group(1)}{org_match.group(1)}{cell_number}"
            return stage, cell_id
    
    # Try EX pattern: Org[a-z]EX[stage][cell]
    ex_match = re.search(r'Org[a-z]EX(\d{2})(\d{4})', filename)
    if ex_match:
        stage = f"EX{ex_match.group(1)}"  # EX05, EX15
        cell_number = ex_match.group(2)
        prefix_match = re.search(r'(Org[a-z]EX)', filename)
        if prefix_match:
            cell_id = f"{prefix_match.group(1)}{ex_match.group(1)}{cell_number}"
            return stage, cell_id
    
    return None, None

def organize_conversion_by_stage(input_dir, output_base_dir, chrom_sizes, resolution=100000, max_files_per_stage=None):
    """
    Convert pairs files organized by developmental stage
    """
    
    # Find all pairs files
    pairs_files = glob.glob(os.path.join(input_dir, "*.pairs.gz"))
    
    if not pairs_files:
        print(f"No .pairs.gz files found in {input_dir}")
        return
    
    # Group files by stage
    stage_files = {}
    failed_extractions = []
    
    for pairs_file in pairs_files:
        stage, cell_id = extract_stage_and_cell_id(pairs_file)
        if stage and cell_id:
            if stage not in stage_files:
                stage_files[stage] = []
            stage_files[stage].append((pairs_file, cell_id))
        else:
            failed_extractions.append(pairs_file)
    
    print(f"Found {len(pairs_files)} total files")
    print(f"Successfully extracted stage info for {len(pairs_files) - len(failed_extractions)} files")
    print(f"Failed extractions: {len(failed_extractions)}")
    
    if failed_extractions:
        print("Files with failed stage extraction:")
        for f in failed_extractions[:5]:  # Show first 5
            print(f"  {os.path.basename(f)}")
        if len(failed_extractions) > 5:
            print(f"  ... and {len(failed_extractions) - 5} more")
    
    print(f"\\nStages found:")
    for stage in sorted(stage_files.keys()):
        count = len(stage_files[stage])
        print(f"  {stage}: {count} cells")
    
    # Process each stage
    for stage in sorted(stage_files.keys()):
        files_for_stage = stage_files[stage]
        
        if max_files_per_stage:
            files_for_stage = files_for_stage[:max_files_per_stage]
            print(f"\\nProcessing first {max_files_per_stage} files for {stage}")
        else:
            print(f"\\nProcessing all {len(files_for_stage)} files for {stage}")
        
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
                else:
                    print(f"    ✗ Failed to convert")
                    print(f"    Error: {result.stderr}")
            except Exception as e:
                print(f"    ✗ Exception: {e}")
        
        print(f"  Stage {stage} completed\\n")
    
    print("=== Conversion Summary ===")
    for stage in sorted(stage_files.keys()):
        stage_dir = os.path.join(output_base_dir, stage)
        if os.path.exists(stage_dir):
            cell_dirs = [d for d in os.listdir(stage_dir) if os.path.isdir(os.path.join(stage_dir, d))]
            print(f"{stage}: {len(cell_dirs)} cells converted")

def main():
    parser = argparse.ArgumentParser(description='Stage-aware conversion of pairs files')
    parser.add_argument('--input_dir', required=True, help='Directory containing .pairs.gz files')
    parser.add_argument('--output_dir', default='converted_matrices_by_stage', help='Output base directory')
    parser.add_argument('--chrom_sizes', required=True, help='Chromosome sizes file')
    parser.add_argument('--resolution', type=int, default=100000, help='Resolution in bp')
    parser.add_argument('--max_per_stage', type=int, help='Maximum files per stage (for testing)')
    
    args = parser.parse_args()
    
    organize_conversion_by_stage(
        args.input_dir,
        args.output_dir, 
        args.chrom_sizes,
        args.resolution,
        args.max_per_stage
    )

if __name__ == "__main__":
    main()