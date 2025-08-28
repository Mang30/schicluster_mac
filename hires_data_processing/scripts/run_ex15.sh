#!/bin/bash
# Processing script for developmental stage EX15

STAGE="EX15"
PROJECT_DIR="/Volumes/SumSung500/CSU/0_HiRES/hires_data_processing"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

echo "Starting processing for stage $STAGE at $(date)"

# Activate environment and run processing
eval "$(micromamba shell hook --shell=bash)"
micromamba activate schicluster

cd "$PROJECT_DIR"

# Run processing for this specific stage
python3 scripts/process_hic_by_stage.py --specific_stage "$STAGE" 2>&1 | tee "$LOG_DIR/processing_$STAGE.log"

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
