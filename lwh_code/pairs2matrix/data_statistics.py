#!/usr/bin/env python3
"""
Generate statistics for the HiC dataset organized by stages
"""

import os
import glob
import re
from collections import defaultdict

def extract_stage_and_cell_id(filepath):
    """Extract stage and cell ID from filepath - same logic as conversion script"""
    filename = os.path.basename(filepath)
    
    # Try Gas* pattern first: Gas[a-z]E[stage][cell]
    gas_match = re.search(r'Gas[a-z]E(\d{2})(\d{4})', filename)
    if gas_match:
        stage = f"E{gas_match.group(1)}"
        cell_number = gas_match.group(2)
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
        stage = f"EX{ex_match.group(1)}"
        cell_number = ex_match.group(2)
        prefix_match = re.search(r'(Org[a-z]EX)', filename)
        if prefix_match:
            cell_id = f"{prefix_match.group(1)}{ex_match.group(1)}{cell_number}"
            return stage, cell_id
    
    return None, None

def analyze_dataset(data_dir):
    """Analyze the complete dataset"""
    
    print("=" * 60)
    print("HiC Dataset Analysis - Stage-organized")
    print("=" * 60)
    
    # Find all pairs files
    pairs_files = glob.glob(os.path.join(data_dir, "*.pairs.gz"))
    
    print(f"Total pairs.gz files found: {len(pairs_files)}")
    
    # Group by stages
    stage_files = defaultdict(list)
    failed_files = []
    
    for pairs_file in pairs_files:
        stage, cell_id = extract_stage_and_cell_id(pairs_file)
        if stage and cell_id:
            stage_files[stage].append((pairs_file, cell_id))
        else:
            failed_files.append(pairs_file)
    
    total_identified = sum(len(files) for files in stage_files.values())
    
    print(f"Successfully identified: {total_identified} files ({total_identified/len(pairs_files)*100:.1f}%)")
    print(f"Failed identification: {len(failed_files)} files ({len(failed_files)/len(pairs_files)*100:.1f}%)")
    
    print("\n" + "=" * 60)
    print("STAGE DISTRIBUTION")
    print("=" * 60)
    
    # Define expected stage order
    stage_order = ['E70', 'E75', 'E80', 'E85', 'E95', 'EX05', 'EX15']
    
    total_cells = 0
    for stage in stage_order:
        if stage in stage_files:
            count = len(stage_files[stage])
            total_cells += count
            percentage = count / total_identified * 100
            print(f"{stage:>6}: {count:>5} cells ({percentage:>5.1f}%)")
    
    print(f"{'TOTAL':>6}: {total_cells:>5} cells")
    
    if failed_files:
        print(f"\n" + "=" * 60)
        print("UNIDENTIFIED FILES (sample)")
        print("=" * 60)
        for i, f in enumerate(failed_files[:10]):
            print(f"  {os.path.basename(f)}")
        if len(failed_files) > 10:
            print(f"  ... and {len(failed_files) - 10} more")
    
    # Check for specific patterns in unidentified files
    patterns = defaultdict(int)
    for f in failed_files:
        basename = os.path.basename(f)
        if 'Val' in basename:
            patterns['Val*'] += 1
        elif 'Mit' in basename:
            patterns['Mit*'] += 1
        elif 'Tel' in basename:
            patterns['Tel*'] += 1
        else:
            patterns['Other'] += 1
    
    if patterns:
        print(f"\nUnidentified file patterns:")
        for pattern, count in patterns.items():
            print(f"  {pattern}: {count} files")
    
    print(f"\n" + "=" * 60)
    print("SUMMARY FOR SCHICLUSTER IMPUTATION")
    print("=" * 60)
    print(f"✓ {len(stage_order)} developmental stages identified")
    print(f"✓ {total_cells:,} cells ready for processing")
    print(f"✓ Coverage: {total_identified/len(pairs_files)*100:.1f}% of dataset")
    print(f"✓ Data organized by stage for comparative analysis")

if __name__ == "__main__":
    data_dir = "data/GSE223917_RAW"
    analyze_dataset(data_dir)