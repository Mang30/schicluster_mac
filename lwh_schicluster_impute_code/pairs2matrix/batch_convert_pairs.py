#!/usr/bin/env python3
"""
Batch convert multiple pairs.gz files to schicluster sparse matrix format
"""

import os
import glob
import subprocess
import argparse
from pathlib import Path

def extract_cell_id(filepath):
    """Extract cell ID from filepath"""
    filename = os.path.basename(filepath)
    # Pattern: GSM6998595_GasaE751001.pairs.gz -> GasaE751001
    if '_' in filename:
        return filename.split('_')[1].split('.')[0]
    else:
        return filename.split('.')[0]

def batch_convert(input_dir, output_dir, chrom_sizes, resolution=100000, max_files=None):
    """
    Batch convert all pairs.gz files in input directory
    """
    
    # Find all pairs.gz files
    pairs_files = glob.glob(os.path.join(input_dir, "*.pairs.gz"))
    
    if not pairs_files:
        print(f"No .pairs.gz files found in {input_dir}")
        return
    
    print(f"Found {len(pairs_files)} pairs files")
    
    if max_files:
        pairs_files = pairs_files[:max_files]
        print(f"Processing first {max_files} files")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert each file
    for i, pairs_file in enumerate(pairs_files):
        cell_id = extract_cell_id(pairs_file)
        cell_output_dir = os.path.join(output_dir, cell_id)
        
        print(f"\n[{i+1}/{len(pairs_files)}] Converting {cell_id}...")
        print(f"Input: {pairs_file}")
        print(f"Output: {cell_output_dir}")
        
        # Run conversion script
        cmd = [
            'python', 'lwh_code/convert_pairs_to_matrix.py',
            '--pairs_file', pairs_file,
            '--output_dir', cell_output_dir,
            '--cell_id', cell_id,
            '--resolution', str(resolution),
            '--chrom_sizes', chrom_sizes
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✓ Successfully converted {cell_id}")
            else:
                print(f"✗ Failed to convert {cell_id}")
                print(f"Error: {result.stderr}")
        except Exception as e:
            print(f"✗ Exception converting {cell_id}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Batch convert pairs files to schicluster format')
    parser.add_argument('--input_dir', required=True, help='Directory containing .pairs.gz files')
    parser.add_argument('--output_dir', default='converted_matrices', help='Output directory')
    parser.add_argument('--chrom_sizes', required=True, help='Chromosome sizes file')
    parser.add_argument('--resolution', type=int, default=100000, help='Resolution in bp')
    parser.add_argument('--max_files', type=int, help='Maximum number of files to process (for testing)')
    
    args = parser.parse_args()
    
    batch_convert(
        args.input_dir,
        args.output_dir, 
        args.chrom_sizes,
        args.resolution,
        args.max_files
    )

if __name__ == "__main__":
    main()