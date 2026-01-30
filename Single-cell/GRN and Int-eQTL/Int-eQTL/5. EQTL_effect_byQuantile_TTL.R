#!/usr/bin/env Rscript

#Author: FQ
#Modified: Thong Luong
#Date: December 5th

args = commandArgs(trailingOnly=TRUE)
print(args)

eqtl_index=as.numeric(args[1]) #Gene index

cell_type = "epithelial" # major cell type

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/')

#eqtl_list <- read.table("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/eQTL_list/eQTL_list_symbol.tsv", header=T)
#eqtl_list$eqtl_ID <- paste0(eqtl_list$phenotype_id, "_", eqtl_list$variant_id)

## Prepare sig eQTL from interaction test, then generate heatmap for each Quantiles_pseudotime
res <- read.csv(file=paste0("Q_pseudotime_output/NBME_top_snps.csv"))
res_sig_FDR <- res[res$Inter_Pval_FDR < 0.05, ]
eqtl_list <- res_sig_FDR[, c(1, 3, 2)]
colnames(eqtl_list) <- c("phenotype_id", "variant_id", "symbol")
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
Meta_branch2 <- Meta[Meta$cell_types %in% c("AT2", "Alveolar Transitional Cells", "AT1"), ]
meta_data1 <- merge(Meta_branch2, age_gpc_cov, by="Sample", all.x=T)

idx <- match(rownames(Meta_branch2), colnames(Expre_matrix))
Expre_matrix_branch2=Expre_matrix[, idx]

idx <- match(rownames(Meta_branch2), rownames(Pse_time))
Pse_time_branch2=Pse_time[idx, ]
Pse_time_branch2$Q_time_B2 <- cut(Pse_time_branch2$pseudotime, breaks = quantile(Pse_time_branch2$pseudotime, probs = c(0:6)/6), include.lowest = T, labels = FALSE)

idx <- match(rownames(Meta_branch2), rownames(ePC_trans))
expPC=ePC_trans[idx, ]

SNP_bed2 <- as.data.frame(readRDS('top_eqtl_snps.rds'))
SNP_bed2_uni <- SNP_bed2[, !duplicated(colnames(SNP_bed2))]
idx <- match(eqtl_list$variant_id, colnames(SNP_bed2_uni))
SNP_bed2 <- SNP_bed2_uni[, idx]

all(eqtl_list$variant_id==colnames(SNP_bed2))
colnames(SNP_bed2) <- eqtl_list$eqtl_ID

SNP_bed2$Sample <- rownames(SNP_bed2)
meta_data2 <- merge(meta_data1, SNP_bed2, by="Sample", all.x=T)

for (j in (10*(eqtl_index-1)+1):min(nrow(eqtl_list),(10*eqtl_index))) {
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

	data$IND <- as.factor(data$IND)

	library(lme4)
    library(glmmTMB)
    for (Q_i in 1:6){
    	data1 <- data[data$Q_pseudotime==Q_i, ]
		full_model <- tryCatch(
                {
                lme4::glmer.nb(formula = E ~ G + (1 | IND) + 
                AGE + nUMI + MT + gPC1 + gPC2 + gPC3 + 
                expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + expPC6, 
                data = data1, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
                },
                error=function(e){
                  message("glmer.nb failed, switching to glmmTMB...")
                  glmmTMB::glmmTMB(E ~ G + (1 | IND) + 
                AGE + nUMI + MT + gPC1 + gPC2 + gPC3 + 
                expPC1 + expPC2 + expPC3 + expPC4 + expPC5 + expPC6, family = nbinom2, data = data1) 
                }
           )
        out = summary(full_model)$coeff
        if (length(out)==3){
          out <- out$cond
        }
		colnames(out) <- c("Estimate","Std.Error","zvalue","pval")

	    message('Save results')
	    outprefix = 'Q_pseudotime_output/EQTL_effect_byQuantile_Finite_time/'
	    out1 = data.frame(term = row.names(out), gene = gene_name, snp = SNP_name, Q_pseudotime=Q_i, out) 

        out1 <- out1[, c("term", "gene", "snp", "Q_pseudotime", "Estimate","Std.Error","zvalue","pval")]
	    write.csv(out1, paste0(outprefix, paste('result', cell_type, gene_name, SNP_name, Q_i, sep = '_'), '.csv'), quote = F)
	}
}

