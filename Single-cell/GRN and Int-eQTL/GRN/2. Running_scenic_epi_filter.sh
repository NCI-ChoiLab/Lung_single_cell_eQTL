#!/bin/bash
#sbatch -c30 --mem=200g --gres=lscratch:50 --time=7-00:00:00 scenic_epi_filter.sh
#July 8th

source /data/Choi_lung/TTL/mambaforge/etc/profile.d/conda.sh
conda activate pyscenic3

pyscenic grn databases/Epi_filtered.loom resources/allTFs_hg38.txt -o tmp/adj_epi_filtered.csv --num_workers 30
pyscenic ctx tmp/adj_epi_filtered.csv databases/*.feather --annotations_fname resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname databases/Epi_filtered.loom --output tmp/reg_epi_filtered.csv --mask_dropouts --num_workers 4
pyscenic aucell databases/Epi_filtered.loom tmp/reg_epi_filtered.csv --output tmp/Epi_filtered_aucell.loom --num_workers 4


