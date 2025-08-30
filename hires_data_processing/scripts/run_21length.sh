#!/bin/bash
# Processing script for 21-length chromosome contact table

STAGE="21length"
# Use relative paths from the script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

echo "Starting processing for stage $STAGE at $(date)"

# Activate environment and run processing
eval "$(micromamba shell hook --shell=bash)"
micromamba activate 3_schicluster_python38

cd "$PROJECT_DIR"

# Run processing for this specific stage with 21-chromosome configuration
python3 scripts/process_hic_by_stage.py --specific_stage "$STAGE" \
    --chrom_size_file "mm10_chrom_sizes_with_chrY.txt" 2>&1 | tee "$LOG_DIR/processing_$STAGE.log"

EXIT_CODE=${PIPESTATUS[0]}

if [ $EXIT_CODE -eq 0 ]; then
    echo "Stage $STAGE completed successfully at $(date)"
    touch "$LOG_DIR/completed_$STAGE.flag"
else
    echo "Stage $STAGE failed with exit code $EXIT_CODE at $(date)"
    touch "$LOG_DIR/failed_$STAGE.flag"
fi

micromamba deactivate
exit $EXIT_CODE