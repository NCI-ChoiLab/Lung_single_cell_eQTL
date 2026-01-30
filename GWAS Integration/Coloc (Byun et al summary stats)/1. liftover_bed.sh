#!/bin/bash
#Jan 16h 2025
#sbatch --partition=quick  -c4 --mem=50g --time=03:00:00 liftover_bed.sh

module load ucsc

liftOver hg19_total.bed hg19ToHg38.over.chain.gz hg38_total.bed total_unlifted.bed
#liftOver hg19_luad.bed hg19ToHg38.over.chain.gz hg38_luad.bed luad_unlifted.bed
#liftOver hg19_lusc.bed hg19ToHg38.over.chain.gz hg38_lusc.bed lusc_unlifted.bed
#liftOver hg19_scc.bed hg19ToHg38.over.chain.gz hg38_scc.bed scc_unlifted.bed

