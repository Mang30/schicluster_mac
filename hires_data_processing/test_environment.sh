#!/bin/bash
# Test script to verify environment setup before running full processing

echo "Testing schicluster environment and dependencies..."

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_test() {
    echo -e "${GREEN}[TEST]${NC} $1"
}

print_pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

print_fail() {
    echo -e "${RED}[FAIL]${NC} $1"
}

print_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

FAILED_TESTS=0

# Test 1: Check micromamba/conda
print_test "Checking micromamba availability..."
if command -v micromamba &> /dev/null; then
    print_pass "micromamba found"
    CONDA_CMD="micromamba"
elif command -v conda &> /dev/null; then
    print_pass "conda found (will use instead of micromamba)"
    CONDA_CMD="conda"
else
    print_fail "Neither micromamba nor conda found"
    ((FAILED_TESTS++))
    CONDA_CMD=""
fi

# Test 2: Check if schicluster environment exists
if [ -n "$CONDA_CMD" ]; then
    print_test "Checking schicluster environment..."
    if $CONDA_CMD env list | grep -q "schicluster"; then
        print_pass "schicluster environment exists"
        
        # Test 3: Activate environment and check hicluster
        print_test "Testing hicluster command in schicluster environment..."
        eval "$($CONDA_CMD shell hook --shell=bash)"
        $CONDA_CMD activate schicluster
        
        if command -v hicluster &> /dev/null; then
            print_pass "hicluster command available"
            
            # Test version
            print_test "Checking hicluster version..."
            VERSION=$(hicluster --version 2>&1 || echo "version check failed")
            echo "  Version: $VERSION"
            
        else
            print_fail "hicluster command not found in schicluster environment"
            print_warn "Install with: $CONDA_CMD install -c bioconda schicluster"
            ((FAILED_TESTS++))
        fi
        
        $CONDA_CMD deactivate
        
    else
        print_fail "schicluster environment not found"
        print_warn "Create with: $CONDA_CMD create -n schicluster -c bioconda schicluster"
        ((FAILED_TESTS++))
    fi
fi

# Test 4: Check file paths
print_test "Checking required file paths..."

CHROM_FILE="/Volumes/SumSung500/CSU/0_HiRES/mm10_chrom_sizes_with_chrY.txt"
if [ -f "$CHROM_FILE" ]; then
    print_pass "Chromosome sizes file found"
    echo "  File: $CHROM_FILE"
    echo "  Size: $(wc -l < "$CHROM_FILE") chromosomes"
else
    print_fail "Chromosome sizes file not found: $CHROM_FILE"
    ((FAILED_TESTS++))
fi

DATA_DIR="/Volumes/SumSung500/CSU/0_HiRES/data/hires/GSE223917_RAW_by_real_stages"
if [ -d "$DATA_DIR" ]; then
    print_pass "Data directory found"
    echo "  Directory: $DATA_DIR"
    STAGE_COUNT=$(find "$DATA_DIR" -maxdepth 1 -type d -name "E*" | wc -l)
    echo "  Stages found: $STAGE_COUNT"
else
    print_fail "Data directory not found: $DATA_DIR"
    ((FAILED_TESTS++))
fi

# Test 5: Check contact tables
print_test "Checking contact tables..."
CONTACT_DIR="/Volumes/SumSung500/CSU/0_HiRES/hires_data_processing/contact_tables"
if [ -d "$CONTACT_DIR" ]; then
    TABLE_COUNT=$(ls "$CONTACT_DIR"/contact_table_*.tsv 2>/dev/null | wc -l)
    if [ "$TABLE_COUNT" -gt 0 ]; then
        print_pass "Contact tables found: $TABLE_COUNT files"
        
        # Check a sample table
        SAMPLE_TABLE=$(ls "$CONTACT_DIR"/contact_table_*.tsv | head -1)
        if [ -f "$SAMPLE_TABLE" ]; then
            CELL_COUNT=$(wc -l < "$SAMPLE_TABLE")
            print_pass "Sample table check: $CELL_COUNT cells in $(basename "$SAMPLE_TABLE")"
        fi
    else
        print_fail "No contact tables found in $CONTACT_DIR"
        print_warn "Run: python3 scripts/generate_contact_tables.py"
        ((FAILED_TESTS++))
    fi
else
    print_fail "Contact tables directory not found"
    ((FAILED_TESTS++))
fi

# Test 6: Check disk space
print_test "Checking available disk space..."
OUTPUT_DIR="/Volumes/SumSung500/CSU/0_HiRES/hires_data_processing/outputs"
if [ -d "$OUTPUT_DIR" ]; then
    AVAIL_SPACE=$(df -h "$OUTPUT_DIR" | awk 'NR==2 {print $4}')
    print_pass "Available space in output directory: $AVAIL_SPACE"
    
    # Warn if less than 100GB
    AVAIL_GB=$(df "$OUTPUT_DIR" | awk 'NR==2 {print int($4/1024/1024)}')
    if [ "$AVAIL_GB" -lt 100 ]; then
        print_warn "Available space is less than 100GB. Processing may require significant disk space."
    fi
else
    print_fail "Output directory not found: $OUTPUT_DIR"
    ((FAILED_TESTS++))
fi

# Test 7: Memory check
print_test "Checking system memory..."
if command -v free &> /dev/null; then
    TOTAL_MEM=$(free -g | awk 'NR==2{print $2}')
    AVAIL_MEM=$(free -g | awk 'NR==2{print $7}')
    print_pass "Total memory: ${TOTAL_MEM}GB, Available: ${AVAIL_MEM}GB"
    
    if [ "$AVAIL_MEM" -lt 32 ]; then
        print_warn "Available memory is less than 32GB. Consider reducing batch size."
    fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS memory check
    TOTAL_MEM=$(sysctl -n hw.memsize | awk '{print int($1/1024/1024/1024)}')
    print_pass "Total memory (macOS): ${TOTAL_MEM}GB"
else
    print_warn "Cannot determine memory information on this system"
fi

# Summary
echo ""
echo "==================== TEST SUMMARY ===================="
if [ $FAILED_TESTS -eq 0 ]; then
    print_pass "All tests passed! Environment is ready for processing."
    echo ""
    echo "Ready to run:"
    echo "  ./scripts/run_all_parallel.sh     # Process all stages in parallel"
    echo "  ./scripts/run_e70.sh             # Process single stage"
    echo "  ./scripts/monitor_processing.sh   # Monitor progress"
else
    print_fail "$FAILED_TESTS test(s) failed. Please fix issues before processing."
    exit 1
fi

echo "======================================================"