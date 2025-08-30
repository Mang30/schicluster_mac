#!/usr/bin/env python3
"""
Generate contact table files for each developmental stage
"""

import os
import glob

def generate_contact_table_for_stage(data_dir, stage, output_dir):
    """Generate contact table for a specific developmental stage"""
    stage_path = os.path.join(data_dir, f"stage_{stage}")
    pairs_files = glob.glob(os.path.join(stage_path, "*.pairs.gz"))
    
    if not pairs_files:
        print(f"Warning: No pairs files found in {stage_path}")
        return
    
    contact_table_path = os.path.join(output_dir, f"contact_table_{stage}.tsv")
    
    with open(contact_table_path, 'w') as f:
        for pairs_file in sorted(pairs_files):
            # Resolve the symlink to get the actual file path
            try:
                actual_file_path = os.path.realpath(pairs_file)
                cell_id = os.path.basename(pairs_file).replace('.pairs.gz', '')
                f.write(f"{cell_id}\t{actual_file_path}\n")
            except Exception as e:
                print(f"Warning: Could not resolve symlink for {pairs_file}: {e}")
                # Write the symlink path if we can't resolve it
                cell_id = os.path.basename(pairs_file).replace('.pairs.gz', '')
                f.write(f"{cell_id}\t{pairs_file}\n")
    
    print(f"Generated {contact_table_path} with {len(pairs_files)} cells")
    return contact_table_path

def main():
    # Use relative paths from the project root
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(project_root, "data")
    output_dir = os.path.join(project_root, "hires_data_processing", "contact_tables")
    
    # Developmental stages
    stages = ["E70", "E75", "E80", "E85", "E95", "EX05", "EX15"]
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate contact tables for each stage
    for stage in stages:
        print(f"\nProcessing stage: {stage}")
        contact_table_path = generate_contact_table_for_stage(data_dir, stage, output_dir)
        
        if contact_table_path and os.path.exists(contact_table_path):
            # Print first few lines to verify
            with open(contact_table_path, 'r') as f:
                lines = f.readlines()[:3]
                print(f"Sample entries:")
                for line in lines:
                    print(f"  {line.strip()}")

if __name__ == "__main__":
    main()