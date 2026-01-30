# Single-cell eQTL from Asian never-smokers lung tissues

This repo includes codes to analyze data for our sc-eQTL dataset as part of the manuscript: "Single-cell lung eQTL dataset of Asian never-smokers highlights the roles of alveolar cells in lung cancer etiology". 

To browse gene expression alongside eQTL results across lung cell types, as well as results from our companion isoform-level dataset (using long-read sequencing) please visit (https://appshare.cancer.gov/ISOLUTION/). 

---

## Analysis

[Cell Annotation](https://github.com/ttluong1/Lung_single_cell_eQTL/tree/main/Cell%20Annotation)

Our process for cell type annotation including preprocessing, QC, and subclustering, following the generation of count file using 10x platorm. Processing of our single-cell data was done predominantely using [Seurat](https://satijalab.org/seurat/).

[Single-cell](https://github.com/ttluong1/Lung_single_cell_eQTL/tree/main/Single-cell)

Analyses following cell type annotation including pseudo-bulked single-cell eQTL mapping using [TensorQTL](https://github.com/broadinstitute/tensorqtl), effect size harmonization via [mashr](https://github.com/stephenslab/mashr), and gene regulatory network ([SCENIC](https://github.com/aertslab/SCENIC)) and int-eQTL analysis incorporating pseudotime from [Monocle3](https://cole-trapnell-lab.github.io/monocle3/). 

[GWAS Integration](https://github.com/ttluong1/Lung_single_cell_eQTL/tree/main/GWAS%20Integration)

Nomination of susceptiblity genes from East Asian and multi-ancestry GWAS summary statistics, Shi et al (PMID: 37236969) Byun et al (PMID: 35915169), was done using [coloc](https://chr1swallace.github.io/coloc/) and [TWAS](http://gusevlab.org/projects/fusion/) (Shi et al summmary statistics only).

[Figures](https://github.com/ttluong1/Lung_single_cell_eQTL/tree/main/Figures)

Codes used to generate figures in our manuscript. 

---

## Data
The raw data can be downloaded on GEO accession number (GSE##).
The genotype data used for sc-eQTL mapping can be found on dBGap under accession: phs##.

----

If you use our single-cell eQTL data, please cite the following paper:

Luong T, Yin J, Li B, Shin JH, Sisay E, Mikhail S, Qin F, Anyaso-Samuel S, Kane A, Golden A, Liu J, Zhang Z, Chang YS, Byun J, Han Y, Landi MT, Mancuso N, Banovich N, Rothman N, Amos C, Lan Q, Yu K, Zhang T, Long E, Shi J, Lee JG, Kim EY, and Choi J. Single-cell lung eQTL dataset of Asian never-smokers highlights the roles of alveolar cells in lung cancer etiology. 2026
