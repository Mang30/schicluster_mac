# Hi-C Data Processing with schicluster

This directory contains scripts and configurations for processing single-cell Hi-C data using schicluster (hicluster) across different developmental stages.

## Project Structure

```
hires_data_processing/
├── scripts/                           # Processing scripts
│   ├── generate_contact_tables.py     # Generate contact tables for each stage
│   ├── process_hic_by_stage.py        # Main processing script
│   ├── run_stage_processing.sh        # Wrapper script with environment activation
│   ├── stage_specific_runner.py       # Generate individual stage scripts
│   ├── run_e70.sh, run_e75.sh, etc.  # Individual stage processing scripts
│   ├── run_all_parallel.sh            # Run all stages in parallel
│   └── monitor_processing.sh          # Monitor processing progress
├── contact_tables/                    # Contact table files for each stage
├── outputs/                          # Processing outputs organized by stage
│   ├── E70/impute/100K/
│   ├── E75/impute/100K/
│   └── ... (other stages)
└── logs/                            # Processing logs and status files
```

## Data Overview

- **7 Developmental Stages**: E70, E75, E80, E85, E95, EX05, EX15
- **Total Cells**: ~7,469 cells across all stages
- **File Format**: `.pairs.gz` files (standard pairs format)
- **Chromosome Reference**: `mm10_chrom_sizes_with_chrY.txt` (includes Y chromosome for male cells)

### Cell Counts per Stage
- E70: 557 cells
- E75: 1,870 cells  
- E80: 559 cells
- E85: 1,125 cells
- E95: 1,108 cells
- EX05: 1,128 cells
- EX15: 1,122 cells

## Processing Parameters

- **Resolution**: 100K (100,000 bp)
- **Batch Size**: 1,536 cells
- **CPU per Job**: 30
- **Input Format**: pairs format (--chr1 1 --pos1 2 --chr2 5 --pos2 6)
- **Output Distance**: 500MB
- **Window/Step Size**: 500MB each

## Usage

### Prerequisites

1. **Environment Setup**:
   ```bash
   micromamba activate schicluster
   # or
   conda activate schicluster
   ```

2. **Required Tools**:
   - schicluster/hicluster installed in the environment
   - Python 3.x
   - Standard Unix tools (bash, etc.)

### Processing Options

#### 1. Process Single Stage
```bash
# Process specific stage
cd /Volumes/SumSung500/CSU/0_HiRES/hires_data_processing
./scripts/run_e70.sh

# Or using the main script
./scripts/run_stage_processing.sh --specific_stage E70
```

#### 2. Process All Stages in Parallel
```bash
cd /Volumes/SumSung500/CSU/0_HiRES/hires_data_processing
./scripts/run_all_parallel.sh
```

#### 3. Monitor Processing Progress
```bash
cd /Volumes/SumSung500/CSU/0_HiRES/hires_data_processing
./scripts/monitor_processing.sh
```

#### 4. Custom Parameters
```bash
./scripts/run_stage_processing.sh \\
    --stages E70 E75 \\
    --resolution 50000 \\
    --batch_size 1024 \\
    --cpu_per_job 16
```

## Output Files

After successful processing, each stage will have:

```
outputs/{STAGE}/impute/100K/
├── imputed matrices (HDF5 format)
├── processing logs
└── intermediate files
```

## Monitoring and Logs

- **Individual Logs**: `logs/processing_{STAGE}.log`
- **Status Files**: 
  - `logs/completed_{STAGE}.flag` - successful completion
  - `logs/failed_{STAGE}.flag` - processing failed
- **Master Log**: `logs/schicluster_processing_YYYYMMDD_HHMMSS.log`

## Troubleshooting

### Common Issues

1. **Environment Not Found**:
   ```bash
   micromamba env list  # Check available environments
   micromamba install -c bioconda schicluster  # Install if missing
   ```

2. **Memory Issues**:
   - Reduce `--batch_size` parameter
   - Reduce `--cpu_per_job` parameter

3. **File Path Issues**:
   - Ensure all paths in contact tables are absolute
   - Verify chromosome size file exists

4. **Permission Issues**:
   ```bash
   chmod +x scripts/*.sh  # Make scripts executable
   ```

### Checking Progress

```bash
# Check running processes
ps aux | grep hicluster

# Check disk space in output directories
du -sh outputs/*/

# Monitor log files
tail -f logs/processing_E75.log
```

## Next Steps

After successful imputation:

1. **Concatenate Cells**: Use `hicluster embed-concatcell-chr`
2. **Merge Chromosomes**: Use `hicluster embed-mergechr`
3. **Downstream Analysis**: Clustering, visualization, etc.

## Contact

For issues or questions about the processing pipeline, check the schicluster documentation or processing logs for detailed error messages.