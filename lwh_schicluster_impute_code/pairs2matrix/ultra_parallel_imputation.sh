#!/bin/bash

# Ultra-parallel imputation script with chromosome-level parallelization
# This script maximizes parallelization by processing multiple cells AND chromosomes simultaneously

echo "=== Ultra-parallel schicluster imputation ==="

# Configuration
INPUT_BASE_DIR="converted_matrices_by_stage_20k"
OUTPUT_BASE_DIR="imputed_matrices_by_stage_20k"
CHROM_SIZES="mm10_chrom_sizes.txt"
RESOLUTION=20000

# Ultra-parallel configuration
MAX_PARALLEL_CELLS=${1:-2}      # Number of cells to process in parallel
MAX_PARALLEL_CHROMS=${2:-10}     # Number of chromosomes per cell to process in parallel

# Imputation parameters (keeping quality unchanged)
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
    exit 1
fi

# Create output directory
mkdir -p ${OUTPUT_BASE_DIR}

# Function to process a single chromosome for a cell
process_chromosome() {
    local cell_dir="$1"
    local stage_output_dir="$2"
    local chrom="$3"
    local cell_id=$(basename "$cell_dir")
    
    input_file="${cell_dir}/${cell_id}_${chrom}.txt"
    
    if [ -f "$input_file" ]; then
        echo "[$$] Imputing ${cell_id} ${chrom}..."
        
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
            2>&1 | grep -v "DEBUG" | sed "s/^/[$$] /"
            
        if [ $? -eq 0 ]; then
            echo "[$$] ✓ ${cell_id} ${chrom} completed"
        else
            echo "[$$] ✗ ${cell_id} ${chrom} failed"
        fi
    else
        echo "[$$] ⚠ ${cell_id} ${chrom} input file not found"
    fi
}

# Function to process a single cell (all chromosomes in parallel)
process_cell_parallel() {
    local cell_dir="$1"
    local stage_output_dir="$2"
    local cell_id=$(basename "$cell_dir")
    
    echo "Processing $cell_id with $MAX_PARALLEL_CHROMS parallel chromosomes"
    
    # Create cell output directory
    mkdir -p "${stage_output_dir}/${cell_id}"
    
    # Process chromosomes in parallel
    echo "$CHROMOSOMES" | tr ' ' '\n' | xargs -I {} -P $MAX_PARALLEL_CHROMS bash -c 'process_chromosome "$1" "$2" "$3"' _ "$cell_dir" "$stage_output_dir" {}
    
    echo "✓ Cell ${cell_id} completed"
}

# Export functions
export -f process_chromosome process_cell_parallel
export CHROMOSOMES PAD STD RP RESOLUTION CHROM_SIZES MAX_PARALLEL_CHROMS

# Find all stage directories
STAGE_DIRS=$(find ${INPUT_BASE_DIR} -maxdepth 1 -type d -name "E*")

if [ -z "$STAGE_DIRS" ]; then
    echo "No stage directories found in ${INPUT_BASE_DIR}"
    exit 1
fi

echo "Configuration:"
echo "  Parallel cells: $MAX_PARALLEL_CELLS"
echo "  Parallel chromosomes per cell: $MAX_PARALLEL_CHROMS"
echo "  Total parallel jobs: $((MAX_PARALLEL_CELLS * MAX_PARALLEL_CHROMS))"
echo ""

echo "Found stage directories:"
total_cells=0
for stage_dir in $STAGE_DIRS; do
    stage=$(basename $stage_dir)
    cell_count=$(find $stage_dir -maxdepth 1 -type d -not -path $stage_dir | wc -l)
    total_cells=$((total_cells + cell_count))
    echo "  $stage: $cell_count cells"
done
echo "Total cells to process: $total_cells"
echo ""

# Process each stage
total_start_time=$(date +%s)
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
    
    # Start timing for this stage
    stage_start_time=$(date +%s)
    
    # Process cells in parallel
    echo "$CELL_DIRS" | xargs -I {} -P $MAX_PARALLEL_CELLS bash -c 'process_cell_parallel "$1" "$2"' _ {} $stage_output_dir
    
    # Calculate elapsed time for this stage
    stage_end_time=$(date +%s)
    stage_duration=$((stage_end_time - stage_start_time))
    echo "  Stage ${stage} completed in ${stage_duration}s"
    echo ""
done

total_end_time=$(date +%s)
total_duration=$((total_end_time - total_start_time))

echo "=== Ultra-parallel imputation completed in ${total_duration}s ==="

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