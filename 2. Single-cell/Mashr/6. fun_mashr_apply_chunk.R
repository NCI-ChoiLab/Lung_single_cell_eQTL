#Author: Thong Luong
#Date: Jan 6th 2025

library(rhdf5)
library(mashr)
library(ashr)
library(data.table)
#library(tidyverse)
library(argparse)
source('/data/Choi_lung/TTL/tensor/mashr/utilities.R')

parser <- ArgumentParser()
parser$add_argument("-cn", "--ch_number", type='double', default=1000)  
parser$add_argument("-ch", "--ch_start", type='double', default=50) # dont exceed ch_number - 50
parser$add_argument("-rn", "--random_number", type='double', default=10) # 100, 200 or 300
args <- parser$parse_args()
ch_number <- args$ch_number
ch_start <- args$ch_start
random_number <- as.numeric(args$random_number)
message("chunk start from: ", ch_start)
setwd("/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/FDR") # change the path

# edit the part and use correct fitted model
if (random_number == 200) {
  mashFit <- readRDS("mashr_allp02_200K_fit.rds")
  suffix_rds <- "_200K.rds"
}else if(random_number == 10){
  mashFit <- readRDS("fitted_10K.rds")
  suffix_rds <- "_10K.rds"
}else if(random_number == 50){
  mashFit <- readRDS("mashr_allp02_50K_fit.rds")
  suffix_rds <- "_50K.rds"
}else if(random_number == 300){
  mashFit <- readRDS("mashr_allp02_300K_fit.rds")
  suffix_rds <- "_300K.rds"
}else {
  mashFit <- readRDS("mashr_allp02_fit.rds")
  suffix_rds <- "_100K.rds"
}
message("random number: ", nrow(mashFit$posterior_weights))

message("mixture proportions for different types of covariance matrix:")
print(get_estimated_pi(mashFit))

# error_median <- median(as.matrix(error), na.rm=TRUE)
mashData <- fun_rds_2_mashr("/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/betas.rds","/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/error.rds", max.missing=0) # change the file name here
# mashData <- mash_set_data(as.matrix(betas), 
#                           as.matrix(error), 
#                           zero_Bhat_Shat_reset=error_median)
# m2 <- mash(mashData, g=get_fitted_g(mashFit), fixg=TRUE) # run mash

message("Applying fit mashr model in chunks...")

setwd("/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/FDR/output/")
features <- unique(gsub("\\|.*", "", row.names(mashData$Bhat)))
features <- features[! features == ""]
feature_chunks <- split(features, 
                        cut(seq_along(features), ch_number, labels=FALSE))

# # save_name <- gsub('(.*).\\w+', paste0('\\1chunk_key_all.txt'), 
# #                   "/data/Choi_lung/scLongreads/tensorqtl/isoform_level/output/all.h5")
# capture.output(feature_chunks, file = save_name)
# 
message("chunk start from: ", ch_start)
for(ch in names(feature_chunks)[c((ch_start+1):(ch_start+50))]){
  message("chunk features: ", head(feature_chunks[[ch]]))
  keep <- grepl(paste(feature_chunks[[ch]], collapse="|"), 
                row.names(mashData$Bhat))
  mashData_chunk <- mash_set_data(mashData$Bhat[keep, ],
                                  mashData$Shat[keep, ],
                                  V=mashData$V)
  mashData_chunk <- mash(mashData_chunk, g=get_fitted_g(mashFit), fixg=TRUE) # run mash
  save_name <- gsub('(.*).\\w+', paste0('\\1chunk_allp02', ch, suffix_rds), 
                    "/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/FDR/output/all.h5") # change the path
  saveRDS(mashData_chunk, save_name)
}

message("Done!")

