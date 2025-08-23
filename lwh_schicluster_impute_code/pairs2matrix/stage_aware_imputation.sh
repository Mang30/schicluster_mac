#!/bin/bash

# Stage-aware batch imputation script for schicluster
# This script processes cells organized by developmental stage

# Usage: 
# 1. First activate schicluster environment: micromamba activate schicluster
# 2. Then run: bash lwh_code/stage_aware_imputation.sh

echo "=== Stage-aware schicluster imputation ==="

# Configuration
INPUT_BASE_DIR="converted_matrices_by_stage_20k"
OUTPUT_BASE_DIR="imputed_matrices_by_stage_20k"
CHROM_SIZES="mm10_chrom_sizes.txt"
RESOLUTION=20000

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

# Check if input directory exists
if [ ! -d "$INPUT_BASE_DIR" ]; then
    echo "Error: Input directory $INPUT_BASE_DIR not found"
    echo "Please run stage-aware conversion first"
    exit 1
fi

# Create output directory
mkdir -p ${OUTPUT_BASE_DIR}

# Find all stage directories
STAGE_DIRS=$(find ${INPUT_BASE_DIR} -maxdepth 1 -type d -name "E*")

if [ -z "$STAGE_DIRS" ]; then
    echo "No stage directories found in ${INPUT_BASE_DIR}"
    exit 1
fi

echo "Found stage directories:"
for stage_dir in $STAGE_DIRS; do
    stage=$(basename $stage_dir)
    cell_count=$(find $stage_dir -maxdepth 1 -type d -not -path $stage_dir | wc -l)
    echo "  $stage: $cell_count cells"
done
echo ""

# Process each stage
for stage_dir in $STAGE_DIRS; do
    stage=$(basename $stage_dir)
    echo "=== Processing Stage $stage ==="
    
    # Create stage-specific output directory
    stage_output_dir="${OUTPUT_BASE_DIR}/${stage}"
    mkdir -p $stage_output_dir
    
    # Find all cell directories in this stage
    CELL_DIRS=$(find ${stage_dir} -maxdepth 1 -type d -not -path ${stage_dir})
    
    if [ -z "$CELL_DIRS" ]; then
        echo "  No cell directories found in $stage"
        continue
    fi
    
    cell_count=$(echo "$CELL_DIRS" | wc -l)
    echo "  Processing $cell_count cells in $stage"
    
    # Process each cell in this stage
    cell_num=0
    for cell_dir in $CELL_DIRS; do
        cell_id=$(basename $cell_dir)
        cell_num=$((cell_num + 1))
        
        echo "    [$cell_num/$cell_count] Processing $cell_id"
        
        # Create cell output directory
        mkdir -p ${stage_output_dir}/${cell_id}
        
        # Process each chromosome
        for chrom in $CHROMOSOMES; do
            input_file="${cell_dir}/${cell_id}_${chrom}.txt"
            
            if [ -f "$input_file" ]; then
                echo "      Imputing ${chrom}..."
                
                hicluster impute-cell \
                    --indir ${cell_dir}/ \
                    --outdir ${stage_output_dir}/${cell_id}/ \
                    --cell ${cell_id} \
                    --chrom ${chrom} \
                    --res ${RESOLUTION} \
                    --chrom_file ${CHROM_SIZES} \
                    --pad ${PAD} \
                    --std ${STD} \
                    --rp ${RP} \
                    2>&1 | grep -v "DEBUG"
                    
                if [ $? -eq 0 ]; then
                    echo "        ✓ ${chrom} completed"
                else
                    echo "        ✗ ${chrom} failed"
                fi
            else
                echo "      ⚠ ${chrom} input file not found"
            fi
        done
        
        echo "      Cell ${cell_id} completed"
    done
    
    echo "  Stage ${stage} completed"
    echo ""
done

echo "=== Stage-aware imputation completed ==="

# Summary
echo "=== Imputation Summary ==="
for stage_dir in $STAGE_DIRS; do
    stage=$(basename $stage_dir)
    output_stage_dir="${OUTPUT_BASE_DIR}/${stage}"
    if [ -d "$output_stage_dir" ]; then
        imputed_cells=$(find $output_stage_dir -maxdepth 1 -type d -not -path $output_stage_dir | wc -l)
        echo "$stage: $imputed_cells cells imputed"
    fi
done