#!/bin/bash
# Monitor processing progress

PROJECT_DIR="/Volumes/SumSung500/CSU/0_HiRES/hires_data_processing"
LOG_DIR="$PROJECT_DIR/logs"

while true; do
    clear
    echo "Hi-C Processing Monitor - $(date)"
    echo "=" {1..50} | tr ' ' '='
    echo ""
    
    # Check for completed stages
    COMPLETED=$(ls "$LOG_DIR"/completed_*.flag 2>/dev/null | wc -l)
    FAILED=$(ls "$LOG_DIR"/failed_*.flag 2>/dev/null | wc -l)
    TOTAL=7  # Total number of stages
    RUNNING=$((TOTAL - COMPLETED - FAILED))
    
    echo "Status Overview:"
    echo "  Completed: $COMPLETED"
    echo "  Failed: $FAILED"
    echo "  Running: $RUNNING"
    echo "  Total: $TOTAL"
    echo ""
    
    # Show individual stage status
    STAGES=("E70" "E75" "E80" "E85" "E95" "EX05" "EX15")
    echo "Stage Details:"
    for stage in "${STAGES[@]}"; do
        if [ -f "$LOG_DIR/completed_$stage.flag" ]; then
            echo "  $stage: ✓ Completed"
        elif [ -f "$LOG_DIR/failed_$stage.flag" ]; then
            echo "  $stage: ✗ Failed"
        elif pgrep -f "run_${stage,,}.sh" > /dev/null; then
            echo "  $stage: 🔄 Running"
        else
            echo "  $stage: ⏳ Waiting"
        fi
    done
    
    echo ""
    echo "Press Ctrl+C to exit monitor"
    sleep 5
done
