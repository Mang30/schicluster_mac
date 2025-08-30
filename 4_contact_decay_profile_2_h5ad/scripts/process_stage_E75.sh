#!/bin/bash

# Script to convert contact decay profiles to h5ad format for stage E75

# Set working directory
cd /home/duxuyan/Projects/schicluster_mac

# Create output directory
mkdir -p 4_contact_decay_profile_2_h5ad/outputs

# Run the conversion script
echo "Converting contact decay profiles to h5ad format..."
python 4_contact_decay_profile_2_h5ad/scripts/convert_decay_profile_to_h5ad.py \
    3_create_contact_decay_profile/outputs/stage_E75_100K \
    4_contact_decay_profile_2_h5ad/outputs/stage_E75_100K.h5ad \
    --stage_name "E75"

echo "Conversion completed!"
echo "Output file: 4_contact_decay_profile_2_h5ad/outputs/stage_E75_100K.h5ad"