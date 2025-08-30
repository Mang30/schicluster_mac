#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# 读取h5ad文件
h5ad_file = "/home/duxuyan/Projects/schicluster_mac/4_contact_decay_profile_2_h5ad/outputs_with_metadata/stage_E75_decay_profiles.h5ad"
adata = sc.read(h5ad_file)

print("数据形状:", adata.shape)
print("\n观测值信息:")
print(adata.obs.head())

print("\nCelltype列的唯一值:")
print(adata.obs['Celltype'].unique())

print("\nCelltype列的数据类型:")
print(adata.obs['Celltype'].dtype)

print("\n检查是否有缺失值:")
print(adata.obs['Celltype'].isnull().sum())