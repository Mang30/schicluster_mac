#!/usr/bin/env python3
"""
make_cell_table.py - 依据 Excel 元数据和 stage_files_mapping.csv 生成 cell_table.tsv

输出格式: 两列，无表头，index 为 cell_id（Cellname），第一列为 pairs.gz 绝对路径
若提供第三列需求，可扩展。当前按 schicluster 文档（两列/或一列+index）模式输出。

用法:
  python make_cell_table.py \
    --metadata "/Volumes/SumSung500/CSU/0_HiRES/5_plot_umap/GSE223917_HiRES_emb_metadata.xlsx" \
    --mapping "/Volumes/SumSung500/CSU/0_HiRES/stage_files_mapping.csv" \
    --raw_root "/Volumes/SumSung500/CSU/0_HiRES/data/hires/GSE223917_RAW" \
    --output "/Volumes/SumSung500/CSU/0_HiRES/11pair2h5ad/cell_table.tsv"

说明:
- 优先使用 mapping 的 filepath；若缺失则根据 raw_root 和 filename 拼接。
- 仅包含 Excel 中存在的 Cellname 的条目。
"""
import argparse
import sys
from pathlib import Path
import pandas as pd


def build_cell_table(metadata_path: Path, mapping_path: Path, raw_root: Path) -> pd.DataFrame:
    md = pd.read_excel(metadata_path)
    required = ['Cellname']
    for c in required:
        if c not in md.columns:
            raise ValueError(f"Missing column in metadata: {c}")
    # 只保留需要的列，避免重复并保持顺序
    cells = md[['Cellname']].dropna().drop_duplicates().copy()
    cells.rename(columns={'Cellname': 'cell_id'}, inplace=True)

    # 读取映射
    mp = pd.read_csv(mapping_path)
    if not set(['filename', 'filepath', 'cellname']).issubset(mp.columns):
        raise ValueError("Mapping file must contain columns: filename, filepath, cellname")

    # 用 cellname 连接
    mp2 = mp[['filename', 'filepath', 'cellname']].dropna().drop_duplicates()

    df = cells.merge(mp2, left_on='cell_id', right_on='cellname', how='left')

    # 如果 filepath 缺失，尝试用 raw_root + filename 拼接
    def infer_path(row):
        if isinstance(row.get('filepath'), str) and len(row['filepath']) > 0:
            return row['filepath']
        fn = row.get('filename')
        if isinstance(fn, str) and len(fn) > 0:
            return str(raw_root / fn)
        return None

    df['pairs_path'] = df.apply(infer_path, axis=1)

    # 仅保留存在路径的记录（可选：不检查存在性，保持全量）
    # 这里不做存在性检查，以便快速生成；如需检查可启用
    # mask_exists = df['pairs_path'].apply(lambda p: Path(p).exists())
    # df = df[mask_exists]

    # 组织输出
    out_df = pd.DataFrame({
        'cell_url': df['pairs_path'].values,
    }, index=df['cell_id'].values)

    return out_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--mapping', required=True)
    ap.add_argument('--raw_root', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    metadata_path = Path(args.metadata)
    mapping_path = Path(args.mapping)
    raw_root = Path(args.raw_root)
    output_path = Path(args.output)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    out_df = build_cell_table(metadata_path, mapping_path, raw_root)

    # 写出为 tsv，无表头，index 为 cell_id
    out_df.to_csv(output_path, sep='\t', header=False, index=True)

    # 简要信息
    print(f"Wrote {len(out_df)} rows to {output_path}")
    print(out_df.head().to_string())


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
