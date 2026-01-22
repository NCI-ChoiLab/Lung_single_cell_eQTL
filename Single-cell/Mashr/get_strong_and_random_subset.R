# Get the strong and random QTL signals using significant eQTL as strong

# Thong Luong
# Jan 6th 2025

library(vroom)
library(dplyr)
library(data.table)
library(collapse)
library(rhdf5)

setwd("/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos")

betas = readRDS('betas.rds')
error = readRDS('error.rds')

rownames = rownames(betas)


sig_dirs <- list.dirs("/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output", recursive = FALSE)
sig_dirs <- paste0(sig_dirs, "/top_qtl_results_all_FDR0.05.txt")
sig_dirs <- sig_dirs[file.exists(sig_dirs)]

sig_tests <- vroom(sig_dirs, show_col_types = FALSE, progress = FALSE,
                   col_select=list(feature_id = "phenotype_id",
                                   variant_id = "variant_id")) %>% 
  distinct() %>% fmutate(id = paste(feature_id, variant_id, sep="|"))

message("Saving ", length(intersect(rownames,sig_tests$id)), " strong QTL.")
betas_strong <- betas[intersect(rownames,sig_tests$id), ]
betas_strong <- betas[sig_tests$id, ]
error_strong <- error[intersect(rownames,sig_tests$id), ]
error_strong <- error[sig_tests$id, ]
save_strong = "fdr_top_eqtl.h5"
h5createFile(save_strong)
ncol <- ifelse(ncol(betas_strong) >= 10, 10, ncol(betas_strong))
nrow <- ifelse(nrow(betas_strong) > 1e4, nrow(betas_strong)/100, nrow(betas_strong)/10)
h5createDataset(file = save_strong, dataset = "betas", dims = dim(betas_strong), 
                chunk = c(1000, ncol))
h5createDataset(file = save_strong, dataset = "error", dims = dim(error_strong), 
                chunk = c(1000, ncol))

h5write(as.matrix(betas_strong), save_strong, "betas")
h5write(as.matrix(error_strong), save_strong, "error")
h5write(rownames(betas_strong), save_strong, "rownames")
h5write(colnames(betas_strong), save_strong, "colnames")


random <- sample(1:nrow(betas), 10000)
betas_random <- betas[random, ]
error_random <- error[random, ]

save_random <- "random_all_10K.h5"
h5createFile(save_random)
ncol <- ifelse(ncol(betas_random) >= 10, 10, ncol(betas_random))
nrow <- ifelse(nrow(betas_random) > 1e4, nrow(betas_random)/100, nrow(betas_random)/10)
h5createDataset(file = save_random, dataset = "betas", dims = dim(betas_random),
                chunk = c(1000, ncol))
h5createDataset(file = save_random, dataset = "error", dims = dim(error_random),
                chunk = c(1000, ncol))

h5write(as.matrix(betas_random), save_random, "betas")
h5write(as.matrix(error_random), save_random, "error")
h5write(rownames(betas_random), save_random, "rownames")
h5write(colnames(betas_random), save_random, "colnames")




message("Done!")


