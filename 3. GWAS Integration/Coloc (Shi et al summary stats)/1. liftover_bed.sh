#!/bin/bash
#sbatch --partition=quick  -c4 --mem=50g --time=03:00:00 liftover_bed.sh


module load ucsc

liftOver asn_snps.bed hg19ToHg38.over.chain.gz hg38_asn_snps.bed hg38_asn_snps_unlifted.bed

