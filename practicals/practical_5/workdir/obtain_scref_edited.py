#### Script to create a subset for the single-cell reference dataset from the Whole Mouse Brain dataset from Allen Brain Atlas
## BM 23-09-2024

import os, re
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import random
import numpy as np
import anndata as ad
from collections import Counter

from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# Specify download base path
download_base = Path('./data/abc_download_root')  # Path to where you would like to write the downloaded data.
abc_cache = AbcProjectCache.from_s3_cache(download_base)


files = ['WMB-10Xv2-CTXsp', 'WMB-10Xv2-HPF', 
         'WMB-10Xv2-HY', 'WMB-10Xv2-Isocortex-1',
         'WMB-10Xv2-Isocortex-2', 'WMB-10Xv2-Isocortex-3',
         'WMB-10Xv2-Isocortex-4', 'WMB-10Xv2-MB',
         'WMB-10Xv2-OLF', 'WMB-10Xv2-TH',
         'WMB-10XMulti', # multi
         'WMB-10Xv3-CB', 'WMB-10Xv3-CTXsp', # v3
         'WMB-10Xv3-HPF', 'WMB-10Xv3-HY',
         'WMB-10Xv3-Isocortex-1', 'WMB-10Xv3-Isocortex-2',
         'WMB-10Xv3-MB', 'WMB-10Xv3-MY',
         'WMB-10Xv3-OLF', 'WMB-10Xv3-P',
         'WMB-10Xv3-PAL', 'WMB-10Xv3-PAL',
         'WMB-10Xv3-STR', 'WMB-10Xv3-TH']

#files = ['WMB-10Xv2-CTXsp',  'WMB-10Xv2-HPF', 'WMB-10XMulti','WMB-10Xv3-CB' ]

# Download metadata
print("Downloading metadata...")
path_metadata = abc_cache.get_directory_metadata('WMB-10X')
print("WMB metadata files:\n\t", path_metadata)

# Get cell and gene metadata
print("Reading cell and gene metadata...")
cell = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='cell_metadata').set_index('cell_label')
gene = abc_cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene').set_index('gene_identifier')


# Iterate over files and print the progress
forced = False

for idx,file in enumerate(files):

    print(f"Processing file {idx + 1}/{len(files)}: {file}")
    
    outfile = f'{download_base}/new_expression_matrices/{file}.h5ad'
    if os.path.exists(outfile) and not forced:
        print(f'{outfile} already exists')
        continue
    

    f = file.split("-")
    d = f[0] +"-" + f[1]
    file2 = abc_cache.get_data_path(directory=d, file_name=f"{file}/raw")
    scadata_file = ad.read_h5ad(file2, backed='r')
    
#    scadata_file = ad.read_h5ad(f"{download_base}/expression_matrices/{file}", backed = 'r')

    # Subset metadata to match adatasc indices
    print(f"Subsetting cell metadata to match {scadata_file.shape[0]} cells in the AnnData object...")
    cell_sub = cell[cell.index.isin(scadata_file.obs.index)]

    # Subset adatasc to only include matching cells
    print("Ensuring that the AnnData object contains only matching cells from the metadata...")
    scadata_file = scadata_file[scadata_file.obs.index.isin(cell_sub.index)].to_memory()

    # Reorder cell_sub to match the index order of adatasc
    print("Reordering cell metadata to match the order of cells in the AnnData object...")
    cell_sub = cell_sub.reindex(scadata_file.obs.index)

    scadata_file.obs = cell_sub
    scadata_file.var = gene

    print("Reading spatial object...")
    adata_spatial = ad.read_h5ad('./data/section80_data.h5ad')

    # Find common cluster-alias
    common_clusters = set(adata_spatial.obs['cluster_alias'].unique()).intersection(set(scadata_file.obs['cluster_alias'].unique()))
    print(f"Length of common clusters: {len(common_clusters)}")

    # Filter adatasc to keep only cluster_alias present in common_clusters
    scadata_file = scadata_file[scadata_file.obs['cluster_alias'].isin(common_clusters)]

    print(f"Number of cells in the single-cells after filtering for cluster_alias: {scadata_file.shape[0]}")

    # Save the full AnnData object
    print(f"Saving the concatenated AnnData object with {scadata_file.shape[0]} cells...")
    scadata_file.write_h5ad( outfile) 

    print(f"Process complete for {file}")


# read each file and subsample cells from each cluster equally, taking 1% from each cluster.
alldata = {}
for file in files:
    print(f"Processing file {file}")
    infile = f'{download_base}/new_expression_matrices/{file}.h5ad'
    adata = ad.read_h5ad(infile)
    sc.pp.subsample(adata, fraction=0.01)
    alldata[file] = adata

    


# Concatenate the adata objects
print("Concatenating all datasets...")
adatasc = ad.concat(alldata)

#adatasc = ad.concat([adata_wmb_CTXspv2,adata_wmb_HPFv2,adata_wmb_HYv2,adata_wmb_Isocortex1v2,adata_wmb_Isocortex2v2,
#                     adata_wmb_Isocortex3v2,adata_wmb_Isocortex4v2,adata_wmb_MBv2,adata_wmb_OLFv2,adata_wmb_THv2,
#                     adata_wmb_CBv3, adata_wmb_CTXspv3, adata_wmb_HPFv3,adata_wmb_HYv3,adata_wmb_Isocortex1v3,
#                     adata_wmb_Isocortex2v3, adata_wmb_MBv3, adata_wmb_MYv3, adata_wmb_OLFv3,
#                     adata_wmb_Pv3, adata_wmb_PALv3, adata_wmb_STRv3, adata_wmb_THv3,
#                    adata_wmb_multi])
print(f"Number of cells in the single-cells after filtering for cluster_alias: {adatasc.shape[0]}")


# Save the full AnnData object
print(f"Saving the concatenated AnnData object with {adatasc.shape[0]} cells...")
adatasc.write_h5ad('./data/sc_adata_commonclusters.h5ad')

# Subset 1% of the cells from the concatenated data
#n_cells = adatasc.shape[0]
#print(f"Subsetting 1% of {n_cells} cells from the full dataset...")
#random_indices = np.random.choice(n_cells, size=int(0.01 * n_cells), replace=False)
#adatasc = adatasc[random_indices].copy()


# Subset metadata to match adatasc indices
print(f"Subsetting cell metadata to match {adatasc.shape[0]} cells in the AnnData object...")
cell_sub = cell[cell.index.isin(adatasc.obs.index)]

# Subset adatasc to only include matching cells
print("Ensuring that the AnnData object contains only matching cells from the metadata...")
adatasc = adatasc[adatasc.obs.index.isin(cell_sub.index)]

# Reorder cell_sub to match the index order of adatasc
print("Reordering cell metadata to match the order of cells in the AnnData object...")
cell_sub = cell_sub.reindex(adatasc.obs.index)

adatasc.obs = cell_sub
adatasc.var = gene

# Save the subsetted data
print(f"Saving the 1% subset of the AnnData object with {adatasc.shape[0]} cells...")
adatasc.write_h5ad('./data/sc_adata_1percent.h5ad')


print("Process complete!")

