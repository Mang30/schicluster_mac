#!/usr/bin/env python3
"""
Convert contact decay profile CSV files to h5ad format
"""

import os
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
import argparse
from pathlib import Path

def read_decay_profile(file_path):
    """Read decay profile CSV file and return a DataFrame"""
    df = pd.read_csv(file_path)
    return df

def create_anndata_from_decay_profiles(input_dir, stage_name):
    """Create AnnData object from all decay profiles in a stage directory"""
    
    # Get all sample directories
    sample_dirs = [d for d in Path(input_dir).iterdir() if d.is_dir()]
    
    # Read all decay profiles
    all_data = []
    sample_names = []
    
    for sample_dir in sample_dirs:
        sample_name = sample_dir.name
        csv_file = sample_dir / f"{sample_name}_decay_profile.csv"
        
        if csv_file.exists():
            df = read_decay_profile(csv_file)
            # Use distance as index and contact_frequency as values
            profile_data = df.set_index('distance')['contact_frequency']
            all_data.append(profile_data)
            sample_names.append(sample_name)
        else:
            print(f"Warning: CSV file not found for {sample_name}")
    
    if not all_data:
        raise ValueError("No valid decay profile files found")
    
    # Combine all data into a matrix
    # Align all profiles to the same distance indices
    combined_df = pd.concat(all_data, axis=1)
    combined_df.columns = sample_names
    
    # Fill NaN values with 0
    combined_df = combined_df.fillna(0)
    
    # Create AnnData object
    # Transpose so samples are rows (observations) and distances are columns (variables)
    adata = ad.AnnData(X=combined_df.T.values)
    
    # Set observation names (samples)
    adata.obs_names = sample_names
    
    # Set variable names (distances)
    adata.var_names = [str(d) for d in combined_df.index]
    
    # Add distance information to var
    adata.var['distance'] = combined_df.index.values
    adata.var['distance_kb'] = combined_df.index.values / 1000.0
    
    # Add log-transformed values if they exist in the original data
    if 'log_distance' in df.columns:
        adata.var['log_distance'] = np.log10(combined_df.index.values)
    
    # Add stage information to obs
    adata.obs['stage'] = stage_name
    
    return adata

def main():
    parser = argparse.ArgumentParser(description="Convert contact decay profiles to h5ad format")
    parser.add_argument("input_dir", help="Input directory containing stage subdirectories with decay profiles")
    parser.add_argument("output_file", help="Output h5ad file path")
    parser.add_argument("--stage_name", default="EX05", help="Stage name for metadata")
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Converting decay profiles from {args.input_dir}")
    print(f"Stage name: {args.stage_name}")
    print(f"Output file: {args.output_file}")
    
    # Create AnnData object
    adata = create_anndata_from_decay_profiles(args.input_dir, args.stage_name)
    
    # Save to h5ad
    adata.write(args.output_file)
    
    print(f"Successfully converted {adata.n_obs} samples with {adata.n_vars} distance points")
    print(f"AnnData shape: {adata.shape}")
    print(f"Observations (samples): {list(adata.obs_names[:5])}...")
    print(f"Variables (distances): {list(adata.var_names[:5])}...")

if __name__ == "__main__":
    main()