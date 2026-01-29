
#=======================
# python script utilized to run nsforest to identify 'necessary and sufficient' markers for cell types
#=======================

import pandas as pd
import scanpy as sc
import nsforest as nf
import os
import logging

# configure logging
logging.basicConfig(filename='nsforest.log', level=logging.INFO)

# set the working directory and paths
working_directory = ''
os.chdir(working_directory)
input_file = 'wholegroupjuly.h5ad'
output_folder = working_directory

# check if input file exists
if not os.path.exists(input_file):
    logging.error(f"Input file '{input_file}' not found.")
    raise FileNotFoundError(f"Input file '{input_file}' not found.")

# create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    logging.info(f"Output folder '{output_folder}' created.")

# read the h5ad file using scanpy
adata = sc.read_h5ad(input_file)

# define cluster header
cluster_header = 'cell_types'

# ensure the cluster_header is present in the adata.obs
if cluster_header not in adata.obs:
    logging.error(f"Cluster header '{cluster_header}' not found in adata.obs.")
    raise KeyError(f"Cluster header '{cluster_header}' not found in adata.obs.")

# call NSForest with all genes evaluated and top 10 markers saved
NSForest_results = nf.NSForest(
    adata,
    cluster_header=cluster_header,
    n_trees=100,
    n_genes_eval=adata.shape[1],  # evaluate all genes
    n_binary_genes=10,  # top 10 markers
)

# print results to log
print(NSForest_results)

# save results to CSV files
NSForest_results.to_csv(os.path.join(output_folder, 'wholegroupjulynsforest_results.csv'), index=False)

# extract and save markers
markers = dict(zip(NSForest_results['clusterName'], NSForest_results['NSForest_markers']))

# save supplementary markers CSV
supplementary_markers = pd.DataFrame.from_dict(markers, orient='index').stack().reset_index()
supplementary_markers.columns = ['cluster', 'marker_index', 'marker']
supplementary_markers.to_csv(os.path.join(output_folder, 'wholegroupnsforest_supplementary_markers.csv'), index=False)

