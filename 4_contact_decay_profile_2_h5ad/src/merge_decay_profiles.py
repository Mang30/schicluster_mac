#!/usr/bin/env python3
"""
Merge single-cell HiC decay profiles into h5ad files for each stage.
"""

import os
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
from typing import List, Dict
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def read_decay_profile(file_path: str) -> pd.DataFrame:
    """
    Read a single decay profile CSV file.
    
    Parameters
    ----------
    file_path : str
        Path to the decay profile CSV file
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the decay profile data
    """
    try:
        df = pd.read_csv(file_path)
        return df
    except Exception as e:
        logger.error(f"Error reading {file_path}: {e}")
        return None


def collect_decay_profiles(stage_dir: str) -> Dict[str, pd.DataFrame]:
    """
    Collect all decay profiles for a given stage.
    
    Parameters
    ----------
    stage_dir : str
        Path to the stage directory containing sample subdirectories
        
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping sample names to their decay profile DataFrames
    """
    sample_profiles = {}
    
    # Get all sample directories
    sample_dirs = [d for d in os.listdir(stage_dir) 
                   if os.path.isdir(os.path.join(stage_dir, d))]
    
    logger.info(f"Found {len(sample_dirs)} samples in {stage_dir}")
    
    for sample_dir in sample_dirs:
        sample_path = os.path.join(stage_dir, sample_dir)
        
        # Find the decay profile CSV file
        csv_files = [f for f in os.listdir(sample_path) 
                     if f.endswith('_decay_profile.csv')]
        
        if not csv_files:
            logger.warning(f"No decay profile CSV found in {sample_path}")
            continue
            
        if len(csv_files) > 1:
            logger.warning(f"Multiple decay profile CSVs found in {sample_path}, using the first one")
            
        csv_file = csv_files[0]
        csv_path = os.path.join(sample_path, csv_file)
        
        # Read the decay profile
        profile_df = read_decay_profile(csv_path)
        if profile_df is not None:
            sample_profiles[sample_dir] = profile_df
        else:
            logger.warning(f"Failed to read decay profile for {sample_dir}")
    
    return sample_profiles


def create_anndata_from_profiles(sample_profiles: Dict[str, pd.DataFrame]) -> ad.AnnData:
    """
    Create an AnnData object from sample decay profiles.
    
    Parameters
    ----------
    sample_profiles : Dict[str, pd.DataFrame]
        Dictionary mapping sample names to their decay profile DataFrames
        
    Returns
    -------
    ad.AnnData
        AnnData object with samples as observations and distance bins as variables
    """
    if not sample_profiles:
        raise ValueError("No sample profiles provided")
    
    # Get all unique distances from all profiles
    all_distances = set()
    for df in sample_profiles.values():
        all_distances.update(df['distance'].unique())
    
    # Sort distances
    sorted_distances = sorted(list(all_distances))
    
    # Create expression matrix
    n_samples = len(sample_profiles)
    n_distances = len(sorted_distances)
    
    # Initialize expression matrix
    expr_matrix = np.zeros((n_samples, n_distances))
    
    # Fill expression matrix
    sample_names = list(sample_profiles.keys())
    distance_to_index = {dist: i for i, dist in enumerate(sorted_distances)}
    
    for sample_idx, (sample_name, profile_df) in enumerate(sample_profiles.items()):
        for _, row in profile_df.iterrows():
            distance = row['distance']
            contact_freq = row['contact_frequency']
            
            if distance in distance_to_index:
                distance_idx = distance_to_index[distance]
                expr_matrix[sample_idx, distance_idx] = contact_freq
    
    # Create AnnData object
    adata = ad.AnnData(expr_matrix)
    
    # Add sample names as observation names
    adata.obs_names = sample_names
    
    # Add distance information as variable names and metadata
    adata.var_names = [f"distance_{dist}" for dist in sorted_distances]
    adata.var['distance'] = sorted_distances
    adata.var['distance_kb'] = [dist / 1000 for dist in sorted_distances]
    
    # Add metadata to observations
    adata.obs['sample_id'] = sample_names
    
    return adata


def process_stage(stage_name: str, stage_path: str, output_dir: str):
    """
    Process a single stage: collect profiles and create h5ad file.
    
    Parameters
    ----------
    stage_name : str
        Name of the stage
    stage_path : str
        Path to the stage directory
    output_dir : str
        Directory to save the output h5ad file
    """
    logger.info(f"Processing stage: {stage_name}")
    
    # Collect decay profiles
    sample_profiles = collect_decay_profiles(stage_path)
    
    if not sample_profiles:
        logger.warning(f"No valid decay profiles found for stage {stage_name}")
        return
    
    logger.info(f"Collected {len(sample_profiles)} decay profiles for stage {stage_name}")
    
    # Create AnnData object
    adata = create_anndata_from_profiles(sample_profiles)
    
    # Save to h5ad file
    output_file = os.path.join(output_dir, f"{stage_name}_decay_profiles.h5ad")
    adata.write(output_file)
    
    logger.info(f"Saved {stage_name} decay profiles to {output_file}")
    logger.info(f"AnnData shape: {adata.shape}")


def main():
    parser = argparse.ArgumentParser(description="Merge single-cell HiC decay profiles into h5ad files")
    parser.add_argument("--input_dir", required=True, help="Input directory containing stage subdirectories")
    parser.add_argument("--output_dir", required=True, help="Output directory for h5ad files")
    parser.add_argument("--stage", required=False, help="Specific stage to process (e.g., stage_E75). If not provided, all stages will be processed.")
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.stage:
        # Process only the specified stage
        stage_path = os.path.join(args.input_dir, args.stage)
        if os.path.isdir(stage_path):
            logger.info(f"Processing specific stage: {args.stage}")
            process_stage(args.stage, stage_path, args.output_dir)
        else:
            logger.error(f"Stage directory not found: {stage_path}")
            return
    else:
        # Get all stage directories (original behavior)
        stage_dirs = [d for d in os.listdir(args.input_dir) 
                      if os.path.isdir(os.path.join(args.input_dir, d)) and d.startswith('stage_')]
        
        logger.info(f"Found {len(stage_dirs)} stages to process")
        
        # Process each stage
        for stage_dir in stage_dirs:
            stage_path = os.path.join(args.input_dir, stage_dir)
            process_stage(stage_dir, stage_path, args.output_dir)


if __name__ == "__main__":
    main()