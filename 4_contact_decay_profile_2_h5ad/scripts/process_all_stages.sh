#!/bin/bash
# Script to process all stages

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate 3_schicluster_python38

# Create output directory
mkdir -p /home/duxuyan/Projects/schicluster_mac/4_contact_decay_profile_2_h5ad/outputs

# Run the merge script for all stages
python /home/duxuyan/Projects/schicluster_mac/4_contact_decay_profile_2_h5ad/src/merge_decay_profiles.py \
    --input_dir /home/duxuyan/Projects/schicluster_mac/3_create_contact_decay_profile/outputs \
    --output_dir /home/duxuyan/Projects/schicluster_mac/4_contact_decay_profile_2_h5ad/outputs

echo "Finished processing all stages"