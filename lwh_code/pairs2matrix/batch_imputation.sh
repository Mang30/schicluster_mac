#!/bin/bash

# Batch imputation script for schicluster
# This script processes multiple cells and chromosomes

# Usage: 
# 1. First activate schicluster environment: micromamba activate schicluster
# 2. Then run: bash lwh_code/batch_imputation.sh

echo "=== Batch schicluster imputation ==="

# Configuration
INPUT_BASE_DIR="converted_matrices"
OUTPUT_BASE_DIR="imputed_matrices"
CHROM_SIZES="mm10_chrom_sizes.txt"
RESOLUTION=100000

# Imputation parameters
PAD=1
STD=1
RP=0.5

# Chromosomes to process
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"

# Check if schicluster is available
if ! command -v hicluster &> /dev/null; then
    echo "Error: hicluster command not found. Please activate schicluster environment first:"
    echo "micromamba activate schicluster"
    exit 1
fi

# Create output directory
mkdir -p ${OUTPUT_BASE_DIR}

# Find all cell directories
CELL_DIRS=$(find ${INPUT_BASE_DIR} -maxdepth 1 -type d -not -path ${INPUT_BASE_DIR})

if [ -z "$CELL_DIRS" ]; then
    echo "No cell directories found in ${INPUT_BASE_DIR}"
    exit 1
fi

echo "Found cell directories:"
echo "$CELL_DIRS"
echo ""

# Process each cell
for cell_dir in $CELL_DIRS; do
    cell_id=$(basename $cell_dir)
    echo "Processing cell: $cell_id"
    
    # Create cell output directory
    mkdir -p ${OUTPUT_BASE_DIR}/${cell_id}
    
    # Process each chromosome
    for chrom in $CHROMOSOMES; do
        input_file="${cell_dir}/${cell_id}_${chrom}.txt"
        
        if [ -f "$input_file" ]; then
            echo "  Imputing ${chrom}..."
            
            hicluster impute-cell \
                --indir ${cell_dir}/ \
                --outdir ${OUTPUT_BASE_DIR}/${cell_id}/ \
                --cell ${cell_id} \
                --chrom ${chrom} \
                --res ${RESOLUTION} \
                --chrom_file ${CHROM_SIZES} \
                --pad ${PAD} \
                --std ${STD} \
                --rp ${RP} \
                2>&1 | grep -v "DEBUG"
                
            if [ $? -eq 0 ]; then
                echo "    ✓ ${chrom} completed"
            else
                echo "    ✗ ${chrom} failed"
            fi
        else
            echo "  ⚠ ${chrom} input file not found: ${input_file}"
        fi
    done
    
    echo "  Cell ${cell_id} completed"
    echo ""
done

echo "=== Batch imputation completed ==="