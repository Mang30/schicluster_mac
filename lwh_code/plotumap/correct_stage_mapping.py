import gzip
import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# Define constants
RAW_DATA_DIR = '/home/duxuyan/Projects/0_HiRES/data/GSE223917_RAW'
OUTPUT_DIR = '/home/duxuyan/Projects/0_HiRES/data/new_processed_hic'
METADATA_FILE = '/home/duxuyan/Projects/0_HiRES/data/GSE223917_HiRES_emb_metadata.xlsx'

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

def extract_cell_id(filename):
    '''
    Extract cell ID from filename.
    '''
    # Remove '.pairs.gz' suffix and everything before the last underscore
    base_name = os.path.basename(filename)
    cell_id = base_name.replace('.pairs.gz', '')
    return cell_id

def get_stage_from_filename(filename):
    '''
    Extract stage information from filename based on metadata mapping.
    Corrected stage mapping according to metadata file.
    '''
    # Extract the prefix from the filename (e.g., 'OrgdEX052251' from 'GSM7005022_OrgdEX052251.pairs.gz')
    base_name = os.path.basename(filename)
    prefix = base_name.split('_')[1].split('.')[0]  # Extract the part before the first dot
    
    # Map the prefix to the correct stage according to the metadata file
    # The metadata shows stages: 'E75' 'E70' 'E80' 'E95' 'EX15' 'E85' 'EX05'
    stage_mapping = {
        'Orga': 'E70',    # E70-E75
        'Orgb': 'E75',
        'Orgc': 'E80',
        'Orgd': 'E85',    # EX05
        'Orge': 'E95',    # EX15
        'Orgf': 'EX15',
        'Orgg': 'EX15',
        'Orgh': 'EX15',
        'Vala': 'E95',    # EX15
        'Valb': 'EX15',
        'Valc': 'EX15',
        'Vald': 'EX15',
        'Gasa': 'EX15',
        'Gasb': 'EX15',
        'Gasc': 'EX15'
    }
    
    # Extract first 4 characters to determine the mapping key
    key = prefix[:4]
    return stage_mapping.get(key, 'Unknown')

def read_pairs_file(filepath):
    '''
    Read and parse a .pairs.gz file.
    '''
    data = []
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                # Pairs format typically: readID chr1 pos1 chr2 pos2 strand1 strand2
                # We're interested in positions, assuming they're in the same chromosome for simplicity
                # For inter-chromosomal contacts, we would need a different approach
                if len(parts) >= 7 and parts[1] == parts[3]:  # Same chromosome
                    data.append((int(parts[2]), int(parts[4])))  # pos1, pos2
    return data

def create_contact_matrix(contact_data, resolution=20000, max_bins=10000):
    '''
    Create a sparse contact matrix from contact data.
    '''
    if not contact_data:
        return csr_matrix((max_bins, max_bins))
    
    # Convert to numpy array for easier manipulation
    contacts = np.array(contact_data)
    
    # Convert positions to bins
    bin1 = contacts[:, 0] // resolution
    bin2 = contacts[:, 1] // resolution
    
    # Filter out bins that exceed max_bins
    valid_indices = (bin1 < max_bins) & (bin2 < max_bins)
    bin1 = bin1[valid_indices]
    bin2 = bin2[valid_indices]
    
    # Create symmetric matrix by including both (bin1, bin2) and (bin2, bin1)
    rows = np.concatenate([bin1, bin2])
    cols = np.concatenate([bin2, bin1])
    data = np.ones(len(rows), dtype=int)
    
    # Create sparse matrix with fixed size
    matrix = csr_matrix((data, (rows, cols)), shape=(max_bins, max_bins))
    
    # Remove duplicate entries by summing them
    matrix.sum_duplicates()
    
    return matrix

def process_single_cell(filepath, max_bins=10000):
    '''
    Process a single .pairs.gz file and return a contact matrix.
    '''
    print(f"Processing {filepath}")
    contact_data = read_pairs_file(filepath)
    contact_matrix = create_contact_matrix(contact_data, max_bins=max_bins)
    return contact_matrix

def main():
    '''
    Main function to process all files and create h5ad objects with corrected stage mapping.
    '''
    # Read metadata
    print("Reading metadata...")
    metadata_df = pd.read_excel(METADATA_FILE)
    
    # List all pairs files
    print("Listing pairs files...")
    pairs_files = [f for f in os.listdir(RAW_DATA_DIR) if f.endswith('.pairs.gz')]
    
    # Group files by stage using corrected mapping
    stage_files = {}
    for file in pairs_files:
        stage = get_stage_from_filename(file)
        if stage not in stage_files:
            stage_files[stage] = []
        stage_files[stage].append(file)
    
    # Print stage distribution
    print("\nStage distribution:")
    for stage, files in stage_files.items():
        print(f"{stage}: {len(files)} files")
    
    # Process each stage
    for stage, files in stage_files.items():
        print(f"\nProcessing stage {stage} with {len(files)} files...")
        
        # Initialize lists to store data
        all_matrices = []
        cell_ids = []
        
        # Process each file
        for file in files:
            filepath = os.path.join(RAW_DATA_DIR, file)
            cell_id = extract_cell_id(file)
            cell_ids.append(cell_id)
            
            # Process the file
            matrix = process_single_cell(filepath)
            all_matrices.append(matrix)
        
        # Combine matrices into a single AnnData object
        if all_matrices:
            try:
                # Convert sparse matrices to dense and stack them
                # Using string conversion to avoid the implicit conversion error
                dense_matrices = [m.toarray().astype(str) for m in all_matrices]
                combined_matrix = np.stack(dense_matrices, axis=0)
                
                # Convert back to numeric
                combined_matrix = combined_matrix.astype(float)
                
                # Flatten the last two dimensions to create a 2D matrix (cells x features)
                # Each feature represents a unique bin pair
                flattened_matrix = combined_matrix.reshape(combined_matrix.shape[0], -1)
                
                # Create variable names (BinPair_0, BinPair_1, ...)
                var_names = [f'BinPair_{i}' for i in range(flattened_matrix.shape[1])]
                
                # Create AnnData object
                adata = sc.AnnData(X=flattened_matrix)
                adata.var_names = var_names
                adata.obs_names = cell_ids
                
                # Add metadata
                # Map cell IDs to metadata
                metadata_subset = metadata_df[metadata_df['Cellname'].isin(cell_ids)]
                
                # Create a series with cell IDs as index for easy mapping
                metadata_indexed = metadata_subset.set_index('Cellname')
                
                # Map cell types to observations
                adata.obs['cell_type'] = metadata_indexed.reindex(cell_ids)['Celltype'].values.astype(str)
                adata.obs['stage'] = stage
                adata.obs['batch'] = '0'  # Convert to string
                adata.obs['n_genes'] = str(flattened_matrix.shape[1])  # Convert to string
                
                # Save to h5ad file
                output_file = os.path.join(OUTPUT_DIR, f'{stage}_hic.h5ad')
                adata.write_h5ad(output_file)
                print(f"Saved {output_file}")
                
            except Exception as e:
                print(f"Error processing stage {stage}: {e}")
                continue

if __name__ == '__main__':
    main()