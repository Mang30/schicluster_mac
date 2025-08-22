#!/usr/bin/env python3
"""
Convert pairs.gz files to schicluster sparse matrix format
"""

import gzip
import os
import sys
from collections import defaultdict
import argparse

def parse_chromosome_sizes(chrom_file):
    """Parse chromosome sizes file"""
    chrom_sizes = {}
    with open(chrom_file, 'r') as f:
        for line in f:
            chrom, size = line.strip().split('\t')
            chrom_sizes[chrom] = int(size)
    return chrom_sizes

def convert_pairs_to_matrix(pairs_file, output_dir, cell_id, resolution, chrom_sizes):
    """
    Convert pairs file to sparse matrix format for each chromosome
    
    Format for schicluster: bin1 bin2 count (0-indexed bins)
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Dictionary to store contacts per chromosome
    contacts = defaultdict(lambda: defaultdict(int))
    
    print(f"Processing {pairs_file}...")
    
    # Read pairs file
    with gzip.open(pairs_file, 'rt') if pairs_file.endswith('.gz') else open(pairs_file, 'r') as f:
        for line_num, line in enumerate(f):
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
                
            # pairs format: readname chr1 pos1 chr2 pos2 ...
            chr1, pos1, chr2, pos2 = parts[1], int(parts[2]), parts[3], int(parts[4])
            
            # Skip if chromosomes not in our reference
            if chr1 not in chrom_sizes or chr2 not in chrom_sizes:
                continue
                
            # Only process intra-chromosomal contacts for now
            if chr1 == chr2:
                # Convert positions to bins (0-indexed)
                bin1 = pos1 // resolution
                bin2 = pos2 // resolution
                
                # Make sure bins are within chromosome bounds
                max_bin = chrom_sizes[chr1] // resolution
                if bin1 < max_bin and bin2 < max_bin:
                    # Store contact (always put smaller bin first)
                    if bin1 <= bin2:
                        contacts[chr1][(bin1, bin2)] += 1
                    else:
                        contacts[chr1][(bin2, bin1)] += 1
            
            if line_num % 100000 == 0:
                print(f"Processed {line_num} lines...")
    
    # Write output files for each chromosome
    total_contacts = 0
    for chrom in contacts:
        output_file = os.path.join(output_dir, f"{cell_id}_{chrom}.txt")
        with open(output_file, 'w') as f:
            for (bin1, bin2), count in contacts[chrom].items():
                f.write(f"{bin1}\t{bin2}\t{count}\n")
                total_contacts += count
        
        print(f"Wrote {len(contacts[chrom])} contacts for {chrom} to {output_file}")
    
    print(f"Total contacts processed: {total_contacts}")
    print(f"Chromosomes with contacts: {list(contacts.keys())}")

def main():
    parser = argparse.ArgumentParser(description='Convert pairs file to schicluster matrix format')
    parser.add_argument('--pairs_file', required=True, help='Input pairs.gz file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--cell_id', required=True, help='Cell ID for naming')
    parser.add_argument('--resolution', type=int, default=100000, help='Resolution in bp (default: 100000)')
    parser.add_argument('--chrom_sizes', required=True, help='Chromosome sizes file')
    
    args = parser.parse_args()
    
    # Parse chromosome sizes
    chrom_sizes = parse_chromosome_sizes(args.chrom_sizes)
    print(f"Loaded {len(chrom_sizes)} chromosomes from {args.chrom_sizes}")
    
    # Convert pairs to matrix
    convert_pairs_to_matrix(
        args.pairs_file, 
        args.output_dir, 
        args.cell_id, 
        args.resolution, 
        chrom_sizes
    )

if __name__ == "__main__":
    main()