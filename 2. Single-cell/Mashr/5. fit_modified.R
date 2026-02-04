#Author: Thong Luong
#Nov 18th 2024

library(rhdf5)
library(mashr)
library(ashr)
library(data.table)
library(tidyverse)
library(argparse)
library(matrixcalc)
source('/data/Choi_lung/TTL/tensor/mashr/utilities.R')

setwd('/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/FDR')

strongTmp <- fun_h5_2_mashr('fdr_top_eqtl.h5')
randomTmp <- fun_h5_2_mashr('/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/random_all_10K.h5')
  
strongTmp$Shat <- abs(strongTmp$Shat)
randomTmp$Shat <- abs(randomTmp$Shat)
  
data.temp <- mash_set_data(randomTmp$Bhat, randomTmp$Shat)
Vhat <- estimate_null_correlation_simple(data.temp)
rm(data.temp)
  
data.strong <- mash_set_data(strongTmp$Bhat, strongTmp$Shat, V=Vhat)
data.random <- mash_set_data(randomTmp$Bhat, randomTmp$Shat, V=Vhat)
  
U.pca = cov_pca(data.strong,5)
U.ed = cov_ed(data.strong, U.pca)
  
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)

saveRDS(m,'fitted_10K.rds')

