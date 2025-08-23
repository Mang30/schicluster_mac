#!/bin/bash

# Parallel Stage-aware batch imputation script for schicluster
# This script processes cells organized by developmental stage with parallel processing

# Usage: 
# 1. First activate schicluster environment: micromamba activate schicluster
# 2. Then run: bash lwh_code/parallel_stage_aware_imputation.sh [MAX_PARALLEL_JOBS]

echo "=== Parallel Stage-aware schicluster imputation ==="

# Configuration
INPUT_BASE_DIR="converted_matrices_by_stage_20k"
OUTPUT_BASE_DIR="imputed_matrices_by_stage_20k"
CHROM_SIZES="mm10_chrom_sizes.txt"
RESOLUTION=20000

# Parallel processing configuration
MAX_PARALLEL_JOBS=${1:-4}  # Default to 4 parallel jobs, can be overridden by command line argument

# Imputation parameters (keeping quality unchanged)
PAD=1
STD=1
RP=0.5

# Chromosomes to process
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"

# Check if GNU parallel is available
if ! command -v parallel &> /dev/null; then
    echo "GNU parallel not found. Installing via micromamba..."
    micromamba install -c conda-forge parallel -y
    if ! command -v parallel &> /dev/null; then
        echo "Error: Could not install GNU parallel. Falling back to xargs method."
        USE_PARALLEL=false
    else
        USE_PARALLEL=true
    fi
else
    USE_PARALLEL=true
fi

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

# Function to process a single cell
process_cell() {
    local cell_dir="$1"
    local stage_output_dir="$2"
    local cell_id=$(basename "$cell_dir")
    local stage=$(basename "$(dirname "$cell_dir")")
    
    echo "    Processing $cell_id in stage $stage (PID: $$)"
    
    # Create cell output directory
    mkdir -p "${stage_output_dir}/${cell_id}"
    
    # Process each chromosome for this cell
    for chrom in $CHROMOSOMES; do
        input_file="${cell_dir}/${cell_id}_${chrom}.txt"
        
        if [ -f "$input_file" ]; then
            echo "      [$$] Imputing ${cell_id} ${chrom}..."
            
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
                2>&1 | grep -v "DEBUG" | sed "s/^/      [$$] /"
                
            if [ $? -eq 0 ]; then
                echo "      [$$] ✓ ${cell_id} ${chrom} completed"
            else
                echo "      [$$] ✗ ${cell_id} ${chrom} failed"
            fi
        else
            echo "      [$$] ⚠ ${cell_id} ${chrom} input file not found"
        fi
    done
    
    echo "    [$$] ✓ Cell ${cell_id} completed"
}

# Export the function for parallel processing
export -f process_cell
export CHROMOSOMES PAD STD RP RESOLUTION CHROM_SIZES

# Find all stage directories
STAGE_DIRS=$(find ${INPUT_BASE_DIR} -maxdepth 1 -type d -name "E*")

if [ -z "$STAGE_DIRS" ]; then
    echo "No stage directories found in ${INPUT_BASE_DIR}"
    exit 1
fi

echo "Found stage directories:"
total_cells=0
for stage_dir in $STAGE_DIRS; do
    stage=$(basename $stage_dir)
    cell_count=$(find $stage_dir -maxdepth 1 -type d -not -path $stage_dir | wc -l)
    total_cells=$((total_cells + cell_count))
    echo "  $stage: $cell_count cells"
done
echo "Total cells to process: $total_cells"
echo "Parallel jobs: $MAX_PARALLEL_JOBS"
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
    echo "  Processing $cell_count cells in $stage with $MAX_PARALLEL_JOBS parallel jobs"
    
    # Start timing for this stage
    stage_start_time=$(date +%s)
    
    if [ "$USE_PARALLEL" = true ]; then
        # Use GNU parallel for better job control
        echo "$CELL_DIRS" | parallel -j $MAX_PARALLEL_JOBS --bar process_cell {} $stage_output_dir
    else
        # Fallback to xargs if parallel is not available
        echo "$CELL_DIRS" | xargs -I {} -P $MAX_PARALLEL_JOBS bash -c 'process_cell "$1" "$2"' _ {} $stage_output_dir
    fi
    
    # Calculate elapsed time for this stage
    stage_end_time=$(date +%s)
    stage_duration=$((stage_end_time - stage_start_time))
    echo "  Stage ${stage} completed in ${stage_duration}s"
    echo ""
done

echo "=== Parallel Stage-aware imputation completed ==="

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