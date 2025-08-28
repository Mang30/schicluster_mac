#!/bin/bash
"""
Wrapper script to run schicluster processing with micromamba environment
"""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

print_error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1"
}

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

print_status "Starting schicluster Hi-C data processing"
print_status "Project directory: $PROJECT_DIR"

# Check if micromamba is available
if ! command -v micromamba &> /dev/null; then
    print_error "micromamba not found. Please install micromamba or use conda instead."
    exit 1
fi

# Activate schicluster environment
print_status "Activating schicluster environment..."
eval "$(micromamba shell hook --shell=bash)"
micromamba activate schicluster

if [ $? -ne 0 ]; then
    print_error "Failed to activate schicluster environment"
    print_error "Please make sure the environment exists: micromamba env list"
    exit 1
fi

# Check if hicluster is available
if ! command -v hicluster &> /dev/null; then
    print_error "hicluster command not found in schicluster environment"
    print_error "Please install schicluster in the environment"
    exit 1
fi

print_status "schicluster environment activated successfully"
print_status "hicluster version: $(hicluster --version 2>&1 || echo 'version check failed')"

# Change to project directory
cd "$PROJECT_DIR" || exit 1

# Run processing script with all arguments passed through
print_status "Starting Hi-C data processing..."
print_status "Command: python3 scripts/process_hic_by_stage.py $@"

python3 scripts/process_hic_by_stage.py "$@"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    print_status "Processing completed successfully!"
else
    print_error "Processing failed with exit code $EXIT_CODE"
fi

print_status "Deactivating environment..."
micromamba deactivate

exit $EXIT_CODE