# Contact Decay Profile to h5ad Converter

This directory contains scripts to convert single-cell HiC contact decay profiles to h5ad format for each developmental stage.

## Directory Structure

- `src/`: Contains the main Python script for merging decay profiles
- `scripts/`: Contains bash scripts for processing each stage
- `outputs/`: Will contain the output h5ad files (created after running the scripts)

## Usage

### Process all stages at once:

```bash
./scripts/process_all_stages.sh
```

### Process individual stages:

```bash
./scripts/process_stage_E75.sh
./scripts/process_stage_E80.sh
```

## Output

The output h5ad files will be saved in the `outputs/` directory, with one file per stage:
- `stage_E75_decay_profiles.h5ad`
- `stage_E80_decay_profiles.h5ad`

Each h5ad file contains:
- Observations (rows): Single cells
- Variables (columns): Distance bins
- Data: Contact frequencies at each distance bin for each cell
- Metadata: Sample IDs and distance information

## Requirements

- micromamba environment `3_schicluster_python38`
- Python packages: scanpy, anndata, pandas, numpy