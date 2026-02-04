#!/usr/bin/env Rscript

#Author: FQ
#Modified by TL, December 2nd 2025

library(Seurat)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/')

args = commandArgs(trailingOnly=TRUE)
print(args)

eqtl_index=as.numeric(args[1]) #Gene index

cell_type = "epithelial" # major cell type

eqtl_list <- read.table("eqtl_list_top.tsv", header=T)
eqtl_list$eqtl_ID <- paste0(eqtl_list$phenotype_id, "_", eqtl_list$variant_id)

cds_monocle3 <- get(load("/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/epithelial_cells_pseudotime_root_Y_1063.RData"))
rm(cds)
Pse_time <- as.data.frame(cds_monocle3@principal_graph_aux@listData$UMAP$pseudotime)
Pse_time$cell_id <- row.names(Pse_time)
colnames(Pse_time) <- c("pseudotime", "cell_id")

SCdata <- readRDS("files_for_FQ_KY_2410/Seurat_file/epithelial_cells.rds")
age_gpc_cov <- as.data.frame(t(read.table("files_for_FQ_KY_2410/Covariates/age_gpc_cov.tsv")))
ePC <- (data.table::fread("files_for_FQ_KY_2410/Covariates/ePCs.tsv", sep=","))
ePC_trans <- t(ePC[, -1])
colnames(ePC_trans) <- ePC$V1

## Expression matrix #####################
#############################################
Countdata <- SCdata@assays$RNA@layers$counts
library(Matrix)
Expre_matrix <- as(Countdata, "dgCMatrix")

############# Normalize expression with log2CPM
##################################################
lib_size <- colSums(Expre_matrix)
# Compute CPM
cpm <- t(t(Expre_matrix) / lib_size) * 1e6
# Take log2
logcpm <- log2(cpm + 1)

Expre_matrix <- logcpm
Meta <- SCdata@meta.data
gene_info <- as.data.frame(SCdata@assays$RNA@features@.Data)
gene_names <- row.names(as.data.frame(SCdata@assays$RNA@features@.Data))
gene_annotation <- data.frame(id=gene_names, gene_short_name=gene_names)
rownames(Expre_matrix) <- gene_names
colnames(Expre_matrix) <- rownames(Meta)

#save(Expre_matrix, file="/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/epithelial_logcpm.RData")
#Expre_matrix <- data.table::fread("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/epithelial_logcpm.txt")

library(dplyr)

message('Make dataframe')
# Make data frame of variables for model
age_gpc_cov$Sample <- rownames(age_gpc_cov)
#meta_data <- reference$meta_data

###### Only cell types in branch 2 were utlized in interaction test 
##########################################################################
Meta_branch2 <- Meta[Meta$cell_types %in% c("AT2", "Alveolar Transitional Cells", "AT1"), ]
#Pse_time_finite = Pse_time[is.finite(Pse_time$pseudotime),]
#finite_cells = intersect(rownames(Meta_branch2),rownames(Pse_time_finite))
#Meta_branch2 = Meta_branch2[finite_cells,]
meta_data1 <- merge(Meta_branch2, age_gpc_cov, by="Sample", all.x=T)

idx <- match(rownames(Meta_branch2), colnames(Expre_matrix))
Expre_matrix_branch2=Expre_matrix[, idx]

idx <- match(rownames(Meta_branch2), rownames(Pse_time))
Pse_time_branch2=Pse_time[idx, ]

idx <- match(rownames(Meta_branch2), rownames(ePC_trans))
expPC=ePC_trans[idx, ]

SNP_bed2 <- as.data.frame(readRDS('top_eqtl_snps.rds'))
all(eqtl_list$variant_id==colnames(SNP_bed2))

colnames(SNP_bed2) <- eqtl_list$eqtl_ID

SNP_bed2$Sample <- rownames(SNP_bed2)
meta_data2 <- merge(meta_data1, SNP_bed2, by="Sample", all.x=T)

for (j in (40*(eqtl_index-1)+1):min(nrow(eqtl_list),(40*eqtl_index))) {
    print(paste0(j, "th eqtl"))

	eqtl_name <- colnames(meta_data2)[j+13]
    SNP_name <- eqtl_list$variant_id[eqtl_list$eqtl_ID==eqtl_name]
    gene_name <- eqtl_list$symbol[eqtl_list$eqtl_ID==eqtl_name]

    E = as.numeric(Expre_matrix_branch2[gene_name, ]) %>% round() # round to get integers
	G = meta_data2[, j+13]
	IND = meta_data2$Sample
	AGE = scale(meta_data2$Age)
	nUMI = scale(log(meta_data2$nCount_RNA))
	MT = scale(meta_data2$percent.mt)

	#harmonyPC = t(reference$Z_corr)
    
	data = data.frame(E, G, IND, AGE, nUMI, MT,
        gPC1 = meta_data2$gPC1, gPC2 = meta_data2$gPC2, gPC3 = meta_data2$gPC3, 
        expPC1 = expPC[,1], expPC2 = expPC[,2], expPC3 = expPC[,3], expPC4 = expPC[,4], expPC5 = expPC[,5], 
        expPC6 = expPC[,6], expPC7 = expPC[,7], expPC8 = expPC[,8], expPC9 = expPC[,9], expPC10 = expPC[,10], 
        pseudotime=Pse_time_branch2$pseudotime)

    #### There are some inf values in pseudotime varible, remove these cells
    ##########################################################################################
    data <- data[is.finite(data$pseudotime)==T, ]

    ####### Generate quantiles (1-6) (Continuous variable) based on pseudotime from Monocle 3
    ##########################################################################################
    data$Q_pseudotime <- cut(data$pseudotime, breaks = quantile(data$pseudotime, probs = c(0:6)/6), 
                      include.lowest = T, labels = FALSE)

    #data$Q_pseudotime <- as.factor(data$Q_pseudotime)

	data$IND <- as.factor(data$IND)
	message('Fit Full NBME model')
	ptm <- proc.time() # start the stopwatch!
	library(lme4)

	full_model <- lme4::glmer.nb(formula = E ~ G + (1 | IND) + 
                AGE + nUMI + MT + gPC1 + gPC2 + gPC3 + 
                expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + expPC6 + 
                Q_pseudotime + G:Q_pseudotime, 
                data = data, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
	out = summary(full_model)$coeff
	colnames(out) <- c("Estimate","Std.Error","zvalue","pval")

	message('Save results')
	outprefix = '/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/Q_pseudotime_output/'
	out1 = data.frame(gene = gene_name, snp = SNP_name, term = row.names(out), out) 

  #out2 <- merge(out1, out_null, by="term", all.x=T)
  out2 <- out1[, c("term", "gene", "snp", "Estimate","Std.Error","zvalue","pval")]
	write.csv(out2, paste0(outprefix, paste('result', cell_type, gene_name, SNP_name, sep = '_'), '.csv'), quote = F)
	#saveRDS(data, paste0(outprefix, paste('data', cell_type, gene_name, SNP_name, sep = '_'), '.rds'))
}


