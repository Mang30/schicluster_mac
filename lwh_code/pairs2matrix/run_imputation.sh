#!/bin/bash

# Script to run schicluster imputation
# You need to run this after activating the schicluster environment

# Example usage after activating environment:
# micromamba activate schicluster
# bash run_imputation.sh

echo "Running schicluster imputation..."

# Set paths
CELL_ID="GasaE751001"
CHROM="chr1"
INPUT_DIR="converted_matrices/${CELL_ID}"
OUTPUT_DIR="imputed_matrices"
CHROM_SIZES="mm10_chrom_sizes.txt"
RESOLUTION=100000

# Create output directory
mkdir -p ${OUTPUT_DIR}/${CELL_ID}

# Run imputation for one chromosome as test
echo "Running imputation for ${CELL_ID} chromosome ${CHROM}..."

hicluster impute-cell \
    --indir ${INPUT_DIR}/ \
    --outdir ${OUTPUT_DIR}/${CELL_ID}/ \
    --cell ${CELL_ID} \
    --chrom ${CHROM} \
    --res ${RESOLUTION} \
    --chrom_file ${CHROM_SIZES} \
    --pad 1 \
    --std 1 \
    --rp 0.5

echo "Imputation completed for ${CHROM}"
echo ""
echo "To run imputation for all chromosomes, you can loop through:"
echo "for chrom in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX; do"
echo "    hicluster impute-cell --indir ${INPUT_DIR}/ --outdir ${OUTPUT_DIR}/${CELL_ID}/ --cell ${CELL_ID} --chrom \$chrom --res ${RESOLUTION} --chrom_file ${CHROM_SIZES} --pad 1 --std 1 --rp 0.5"
echo "done"