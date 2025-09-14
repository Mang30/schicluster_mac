# UMAP 绘图（按 celltype 着色，按 stage 分面）

本脚本基于 `31_svd_embedding/output/decomp/total_decomp.npz`（3288×50）做 UMAP 降维，用 `cell_labels.csv` 提供的 `celltype` 和 `stage` 作为分组变量：

- 颜色：celltype（group.by="celltype"）
- 分面：stage（split.by="stage"）

## 依赖

- Python ≥ 3.8
- numpy, pandas, seaborn, matplotlib
- umap-learn（如未安装，会退化为 sklearn TSNE）

## 运行

在仓库根目录或 `31_svd_embedding/plot` 下运行：

```bash
python 31_svd_embedding/plot/plot_umap_by_celltype_split_stage.py \
  --embedding 31_svd_embedding/output/decomp/total_decomp.npz \
  --n_neighbors 15 \
  --min_dist 0.1 \
  --out_prefix 31_svd_embedding/plot/output/umap_celltype_split_stage
```

脚本会：

1. 自动检测并生成 `31_svd_embedding/output/decomp/cell_labels.csv`（若不存在），通过 `code/cell_table.tsv` 与 `stage_files_mapping.csv` 关联。
2. 运行 UMAP，将坐标与标签保存为 CSV。
3. 生成按 stage 分面的 PNG/PDF 图。

输出：

- `31_svd_embedding/plot/output/umap_celltype_split_stage.png`
- `31_svd_embedding/plot/output/umap_celltype_split_stage.pdf`
- `31_svd_embedding/plot/output/umap_celltype_split_stage.csv`

可选：如果存在 `5_plot_umap/color_mapping.json`，会优先使用其中对 celltype 的配色。
