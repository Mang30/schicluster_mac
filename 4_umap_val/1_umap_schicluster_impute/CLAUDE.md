# Hi-C UMAP Visualization Analysis

## Project Overview

This directory contains tools for generating UMAP visualizations from imputed single-cell Hi-C data. The analysis processes 557 cells from the E70 developmental stage stored as .cool format files.

## Data Source

- **Input Directory**: `/0_HiRES/hires_data_processing/outputs/E70/impute/100K/chunk0/`
- **Data Format**: Individual .cool files for each cell
- **Cell Count**: 557 cells
- **Resolution**: 100K bins
- **Cell Metadata**: Available in `cell_table.csv`

## Analysis Pipeline

### 1. Data Loading
- Loads Hi-C contact matrices from .cool files using the cooler library
- Extracts upper triangular features from contact matrices
- Handles missing values and normalizes data

### 2. Feature Processing
- Flattens Hi-C contact matrices to feature vectors
- Applies feature selection for high-dimensional data
- Standardizes features for UMAP input

### 3. UMAP Embedding
Tests multiple parameter configurations:
- **Default**: n_neighbors=15, min_dist=0.1, metric=euclidean
- **More Neighbors**: n_neighbors=30, min_dist=0.1, metric=euclidean  
- **Larger Min Distance**: n_neighbors=15, min_dist=0.5, metric=euclidean
- **Cosine Metric**: n_neighbors=15, min_dist=0.1, metric=cosine

### 4. Visualization
- Individual UMAP plots for each parameter set
- Comparison grid showing all configurations
- Cell labeling (subset for clarity)
- High-resolution output (300 DPI)

## Files

- `hic_umap_analysis.py` - Main analysis script
- `requirements.txt` - Python dependencies
- `run_analysis.sh` - Execution script
- `CLAUDE.md` - This documentation
- `results/` - Output directory (created during analysis)

## Output Structure

```
results/
├── umap_default.png
├── umap_more_neighbors.png  
├── umap_larger_mindist.png
├── umap_cosine_metric.png
├── umap_comparison.png
├── umap_coordinates_default.csv
├── umap_coordinates_more_neighbors.csv
├── umap_coordinates_larger_mindist.csv
└── umap_coordinates_cosine_metric.csv
```

## 使用方法

```bash
# 测试环境
eval "$(micromamba shell hook --shell bash)"
micromamba activate schicluster
python test_environment.py

# 运行完整分析
bash run_analysis.sh

# 或直接运行(需要先激活环境)
micromamba activate schicluster
python hic_umap_analysis.py

# 快速测试(处理少量细胞)
python hic_umap_analysis.py --max_cells 10 --output_dir ./test_results
```

## Dependencies

- cooler: Hi-C data manipulation
- numpy, pandas: Data processing
- scikit-learn: Preprocessing and PCA
- umap-learn: UMAP embedding
- matplotlib, seaborn: Visualization

## Notes

- Analysis can be limited to subset of cells for faster testing
- PCA preprocessing applied automatically for high-dimensional data (>5000 features)
- Results include both visualizations and coordinate data for further analysis
- Error handling implemented for corrupted .cool files
- In this project, please reply in Chinese.
- using micromamba activate schicluster