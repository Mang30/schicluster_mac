#!/bin/bash
# Run all developmental stages in parallel

# Use relative paths from the script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
STAGES=("E70" "E75" "E80" "E85" "E95" "EX05" "EX15")

echo "Starting parallel processing of ${#STAGES[@]} developmental stages"
echo "Stages: ${STAGES[@]}"

# Create logs directory
mkdir -p "$PROJECT_DIR/logs"

# Start all stages in background
PIDS=()
for stage in "${STAGES[@]}"; do
    echo "Starting stage $stage..."
    "$PROJECT_DIR/scripts/run_${stage,,}.sh" &
    PIDS+=($!)
    sleep 2  # Small delay to avoid resource conflicts
done

echo "All stages started. PIDs: ${PIDS[@]}"
echo "Monitor progress with: tail -f $PROJECT_DIR/logs/processing_*.log"

# Wait for all processes to complete
SUCCESS_COUNT=0
TOTAL_COUNT=${#PIDS[@]}

for i in "${!PIDS[@]}"; do
    PID="${PIDS[i]}"
    STAGE="${STAGES[i]}"
    
    echo "Waiting for stage $STAGE (PID: $PID)..."
    wait $PID
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "✓ Stage $STAGE completed successfully"
        ((SUCCESS_COUNT++))
    else
        echo "✗ Stage $STAGE failed (exit code: $EXIT_CODE)"
    fi
done

echo ""
echo "===== SUMMARY ====="
echo "Total stages: $TOTAL_COUNT"
echo "Successful: $SUCCESS_COUNT"
echo "Failed: $((TOTAL_COUNT - SUCCESS_COUNT))"

if [ $SUCCESS_COUNT -eq $TOTAL_COUNT ]; then
    echo "🎉 All stages completed successfully!"
    exit 0
else
    echo "⚠️  Some stages failed. Check individual logs for details."
    exit 1
fi
