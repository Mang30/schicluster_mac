#!/usr/bin/env python3
"""
Generate individual processing scripts for each developmental stage
This allows parallel processing of different stages
"""

import os
from pathlib import Path

def create_stage_script(stage, base_dir):
    """Create a processing script for a specific stage"""
    
    script_content = f"""#!/bin/bash
# Processing script for developmental stage {stage}

STAGE="{stage}"
PROJECT_DIR="{base_dir}"
LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

echo "Starting processing for stage $STAGE at $(date)"

# Activate environment and run processing
eval "$(micromamba shell hook --shell=bash)"
micromamba activate schicluster

cd "$PROJECT_DIR"

# Run processing for this specific stage
python3 scripts/process_hic_by_stage.py --specific_stage "$STAGE" 2>&1 | tee "$LOG_DIR/processing_$STAGE.log"

EXIT_CODE=${{PIPESTATUS[0]}}

if [ $EXIT_CODE -eq 0 ]; then
    echo "Stage $STAGE completed successfully at $(date)"
    touch "$LOG_DIR/completed_$STAGE.flag"
else
    echo "Stage $STAGE failed with exit code $EXIT_CODE at $(date)"
    touch "$LOG_DIR/failed_$STAGE.flag"
fi

micromamba deactivate
exit $EXIT_CODE
"""
    
    script_path = os.path.join(base_dir, "scripts", f"run_{stage.lower()}.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(script_path, 0o755)
    print(f"Created processing script: {script_path}")
    return script_path

def create_parallel_runner(stages, base_dir):
    """Create a script to run all stages in parallel"""
    
    script_content = f"""#!/bin/bash
# Run all developmental stages in parallel

PROJECT_DIR="{base_dir}"
STAGES=({' '.join([f'"{stage}"' for stage in stages])})

echo "Starting parallel processing of ${{#STAGES[@]}} developmental stages"
echo "Stages: ${{STAGES[@]}}"

# Create logs directory
mkdir -p "$PROJECT_DIR/logs"

# Start all stages in background
PIDS=()
for stage in "${{STAGES[@]}}"; do
    echo "Starting stage $stage..."
    "$PROJECT_DIR/scripts/run_${{stage,,}}.sh" &
    PIDS+=($!)
    sleep 2  # Small delay to avoid resource conflicts
done

echo "All stages started. PIDs: ${{PIDS[@]}}"
echo "Monitor progress with: tail -f $PROJECT_DIR/logs/processing_*.log"

# Wait for all processes to complete
SUCCESS_COUNT=0
TOTAL_COUNT=${{#PIDS[@]}}

for i in "${{!PIDS[@]}}"; do
    PID="${{PIDS[i]}}"
    STAGE="${{STAGES[i]}}"
    
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
"""
    
    script_path = os.path.join(base_dir, "scripts", "run_all_parallel.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(script_path, 0o755)
    print(f"Created parallel runner script: {script_path}")
    return script_path

def create_monitor_script(base_dir):
    """Create a monitoring script to check processing status"""
    
    script_content = f"""#!/bin/bash
# Monitor processing progress

PROJECT_DIR="{base_dir}"
LOG_DIR="$PROJECT_DIR/logs"

while true; do
    clear
    echo "Hi-C Processing Monitor - $(date)"
    echo "=" {{1..50}} | tr ' ' '='
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
    for stage in "${{STAGES[@]}}"; do
        if [ -f "$LOG_DIR/completed_$stage.flag" ]; then
            echo "  $stage: ✓ Completed"
        elif [ -f "$LOG_DIR/failed_$stage.flag" ]; then
            echo "  $stage: ✗ Failed"
        elif pgrep -f "run_${{stage,,}}.sh" > /dev/null; then
            echo "  $stage: 🔄 Running"
        else
            echo "  $stage: ⏳ Waiting"
        fi
    done
    
    echo ""
    echo "Press Ctrl+C to exit monitor"
    sleep 5
done
"""
    
    script_path = os.path.join(base_dir, "scripts", "monitor_processing.sh")
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(script_path, 0o755)
    print(f"Created monitor script: {script_path}")
    return script_path

def main():
    base_dir = "/Volumes/SumSung500/CSU/0_HiRES/hires_data_processing"
    stages = ["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"]
    
    print("Creating stage-specific processing scripts...")
    
    # Create individual stage scripts
    for stage in stages:
        create_stage_script(stage, base_dir)
    
    # Create parallel runner
    create_parallel_runner(stages, base_dir)
    
    # Create monitor script
    create_monitor_script(base_dir)
    
    print("\nCreated processing infrastructure:")
    print("1. Individual stage scripts: run_e70.sh, run_e75.sh, etc.")
    print("2. Parallel runner: run_all_parallel.sh")
    print("3. Progress monitor: monitor_processing.sh")
    print("\nUsage examples:")
    print(f"  # Process single stage: {base_dir}/scripts/run_e70.sh")
    print(f"  # Process all stages in parallel: {base_dir}/scripts/run_all_parallel.sh")
    print(f"  # Monitor progress: {base_dir}/scripts/monitor_processing.sh")

if __name__ == "__main__":
    main()