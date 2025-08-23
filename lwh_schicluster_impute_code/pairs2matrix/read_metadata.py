#!/usr/bin/env python3
"""
Read metadata from Excel file to understand stage information
"""

import pandas as pd
import sys

def read_metadata(excel_path):
    """Read metadata from Excel file"""
    try:
        # Read the Excel file
        df = pd.read_excel(excel_path)
        
        print(f"Metadata shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        print("\nFirst 10 rows:")
        print(df.head(10))
        
        # Check if there's a Stage column
        stage_cols = [col for col in df.columns if 'stage' in col.lower() or 'Stage' in col]
        if stage_cols:
            print(f"\nStage-related columns: {stage_cols}")
            for col in stage_cols:
                unique_stages = df[col].unique()
                print(f"\nUnique values in {col}:")
                for stage in sorted(unique_stages):
                    count = (df[col] == stage).sum()
                    print(f"  {stage}: {count} cells")
        else:
            print("\nNo obvious stage column found. All columns:")
            for col in df.columns:
                print(f"  {col}")
                if df[col].dtype == 'object':
                    unique_vals = df[col].unique()
                    if len(unique_vals) <= 20:
                        print(f"    Unique values: {unique_vals}")
                    else:
                        print(f"    {len(unique_vals)} unique values")
                
    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python read_metadata.py <excel_file>")
        sys.exit(1)
    
    read_metadata(sys.argv[1])