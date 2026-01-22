#Author: Thong Luong
#Date: Jan 6th 2025

library(dplyr)
library(tidyr)
library(vroom)
library(collapse)
library(rhdf5)

tensor_path <- "/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/"
dirs <- list.dirs(file.path(tensor_path,"output"), full.names = TRUE, recursive =FALSE)
dirs <- paste0(dirs, "/qtl_results_all.txt")
dirs <- dirs[file.exists(dirs)]

nCelltypes <- length(dirs)

all_qtl <- vroom(dirs, show_col_types = FALSE, id="celltype",
                 col_select=list(celltype="celltype", feature_id = "phenotype_id",
                                 variant_id = "variant_id", 
                                 betas = "slope", error = "slope_se", 
                                 pval_nonimal = "pval_nominal"))

all_qtl <- all_qtl %>%
  # drop duplicates that exist if non-biallelic variants were tested
  funique(cols=c("celltype", "feature_id", "variant_id")) %>% 
  fmutate(celltype = basename(dirname(celltype)),
          id = paste(feature_id, variant_id, sep="|"))


betas <- all_qtl %>% pivot_wider(names_from = celltype, values_from = betas,
                                 id_cols = id) %>%
  tibble::column_to_rownames(var = "id") %>% qDF()
colnames(betas)
dim(betas)
error <- all_qtl %>% pivot_wider(names_from = celltype, values_from = error,
                                 id_cols = id) %>%
  tibble::column_to_rownames(var = "id") %>% qDF()
colnames(error)
betas[is.na(betas)] <- 0
k <- which(is.na(error), arr.ind=TRUE)
error[k] <- rowMeans(error, na.rm=TRUE)[k[, 1]] 



message("Snapshot of merged data...")
n <- ifelse(ncol(betas) > 5, 5, ncol(betas))
message("betas:")
betas[1:5, 1:n]

message("beta standard error:")
error[1:5, 1:n]


### Save output
message("Saving merged results as an hdf5...")
setwd("/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos")
saveRDS(betas, 'betas.rds')
saveRDS(error, 'error.rds')

h5createFile("merge.h5")
ncol <- ifelse(ncol(betas) >= 10, 10, ncol(betas))
nrow <- ifelse(nrow(betas) > 1e4, nrow(betas)/100, nrow(betas)/10)
h5createDataset(file = "merge.h5", dataset = "betas", dims = dim(betas), 
                chunk = c(1000, ncol))
h5createDataset(file = "merge.h5", dataset = "error", dims = dim(error), 
                chunk = c(1000, ncol))

h5write(as.matrix(betas), "merge.h5", "betas")
h5write(as.matrix(error), "merge.h5", "error")
h5write(rownames(betas), "merge.h5", "rownames")
h5write(colnames(betas), "merge.h5", "colnames")

