#!/bin/bash

module load python/3.8 R/4.2

python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Adv_fib.bed.gz Adv_fib --covariates ./Covariates_Sum_Final_chr_pos/Adv_fib/Adv_fib_PF12.txt --mode cis -o ./Output_Sum_Final_chr_pos/Adv_fib
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Alv_fib.bed.gz Alv_fib --covariates ./Covariates_Sum_Final_chr_pos/Alv_fib/Alv_fib_PF12.txt --mode cis -o ./Output_Sum_Final_chr_pos/Alv_fib
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Alv_mph.bed.gz Alv_mph --covariates ./Covariates_Sum_Final_chr_pos/Alv_mph/Alv_mph_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Alv_mph
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Alv_trans.bed.gz Alv_trans --covariates ./Covariates_Sum_Final_chr_pos/Alv_trans/Alv_trans_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Alv_trans
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/AT1.bed.gz AT1 --covariates ./Covariates_Sum_Final_chr_pos/AT1/AT1_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/AT1
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/AT2.bed.gz AT2 --covariates ./Covariates_Sum_Final_chr_pos/AT2/AT2_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/AT2
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Basal.bed.gz Basal --covariates ./Covariates_Sum_Final_chr_pos/Basal/Basal_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Basal
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Bcells.bed.gz Bcells --covariates ./Covariates_Sum_Final_chr_pos/Bcells/Bcells_PF6.txt --mode cis -o ./Output_Sum_Final_chr_pos/Bcells
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/CD4.bed.gz CD4 --covariates ./Covariates_Sum_Final_chr_pos/CD4/CD4_PF14.txt --mode cis -o ./Output_Sum_Final_chr_pos/CD4

python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/CD8.bed.gz CD8 --covariates ./Covariates_Sum_Final_chr_pos/CD8/CD8_PF12.txt --mode cis -o ./Output_Sum_Final_chr_pos/CD8
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Cla_mono.bed.gz Cla_mono --covariates ./Covariates_Sum_Final_chr_pos/Cla_mono/Cla_mono_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Cla_mono
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Club.bed.gz Club --covariates ./Covariates_Sum_Final_chr_pos/Club/Club_PF16.txt --mode cis -o ./Output_Sum_Final_chr_pos/Club
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/DC1.bed.gz DC1 --covariates ./Covariates_Sum_Final_chr_pos/DC1/DC1_PF8.txt --mode cis -o ./Output_Sum_Final_chr_pos/DC1
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/DC2.bed.gz DC2 --covariates ./Covariates_Sum_Final_chr_pos/DC2/DC2_PF16.txt --mode cis -o ./Output_Sum_Final_chr_pos/DC2
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/EC_aero_cap.bed.gz EC_aero_cap --covariates ./Covariates_Sum_Final_chr_pos/EC_aero_cap/EC_aero_cap_PF2.txt --mode cis -o ./Output_Sum_Final_chr_pos/EC_aero_cap
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/EC_art.bed.gz EC_art --covariates ./Covariates_Sum_Final_chr_pos/EC_art/EC_art_PF18.txt --mode cis -o ./Output_Sum_Final_chr_pos/EC_art
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/EC_gen_cap.bed.gz EC_gen_cap --covariates ./Covariates_Sum_Final_chr_pos/EC_gen_cap/EC_gen_cap_PF4.txt --mode cis -o ./Output_Sum_Final_chr_pos/EC_gen_cap

python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/EC_ven_pul.bed.gz EC_ven_pul --covariates ./Covariates_Sum_Final_chr_pos/EC_ven_pul/EC_ven_pul_PF10.txt --mode cis -o ./Output_Sum_Final_chr_pos/EC_ven_pul
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/EC_ven_sys.bed.gz EC_ven_sys --covariates ./Covariates_Sum_Final_chr_pos/EC_ven_sys/EC_ven_sys_PF10.txt --mode cis -o ./Output_Sum_Final_chr_pos/EC_ven_sys
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Goblet.bed.gz Goblet --covariates ./Covariates_Sum_Final_chr_pos/Goblet/Goblet_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Goblet
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Int_mph_peri.bed.gz Int_mph_peri --covariates ./Covariates_Sum_Final_chr_pos/Int_mph_peri/Int_mph_peri_PF8.txt --mode cis -o ./Output_Sum_Final_chr_pos/Int_mph_peri
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Lym_EC_mat.bed.gz Lym_EC_mat --covariates ./Covariates_Sum_Final_chr_pos/Lym_EC_mat/Lym_EC_mat_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Lym_EC_mat
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Lym_EC_pro.bed.gz Lym_EC_pro --covariates ./Covariates_Sum_Final_chr_pos/Lym_EC_pro/Lym_EC_pro_PF8.txt --mode cis -o ./Output_Sum_Final_chr_pos/Lym_EC_pro
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Mast.bed.gz Mast --covariates ./Covariates_Sum_Final_chr_pos/Mast/Mast_PF8.txt --mode cis -o ./Output_Sum_Final_chr_pos/Mast
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Mig_DC.bed.gz Mig_DC --covariates ./Covariates_Sum_Final_chr_pos/Mig_DC/Mig_DC_PF2.txt --mode cis -o ./Output_Sum_Final_chr_pos/Mig_DC

python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Mono_mph.bed.gz Mono_mph --covariates ./Covariates_Sum_Final_chr_pos/Mono_mph/Mono_mph_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Mono_mph
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Multiciliated.bed.gz Multiciliated --covariates ./Covariates_Sum_Final_chr_pos/Multiciliated/Multiciliated_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Multiciliated
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/NK.bed.gz NK --covariates ./Covariates_Sum_Final_chr_pos/NK/NK_PF18.txt --mode cis -o ./Output_Sum_Final_chr_pos/NK
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Noncla_mono.bed.gz Noncla_mono --covariates ./Covariates_Sum_Final_chr_pos/Noncla_mono/Noncla_mono_PF16.txt --mode cis -o ./Output_Sum_Final_chr_pos/Noncla_mono
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Peri_fib.bed.gz Peri_fib --covariates ./Covariates_Sum_Final_chr_pos/Peri_fib/Peri_fib_PF2.txt --mode cis -o ./Output_Sum_Final_chr_pos/Peri_fib
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Sec_trans.bed.gz Sec_trans --covariates ./Covariates_Sum_Final_chr_pos/Sec_trans/Sec_trans_PF20.txt --mode cis -o ./Output_Sum_Final_chr_pos/Sec_trans
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/SM.bed.gz SM --covariates ./Covariates_Sum_Final_chr_pos/SM/SM_PF8.txt --mode cis -o ./Output_Sum_Final_chr_pos/SM
python3 -m tensorqtl ./Genotype/PLINK/chr_pos/final ./Phenotype/Sum_Final/Sub_fib.bed.gz Sub_fib --covariates ./Covariates_Sum_Final_chr_pos/Sub_fib/Sub_fib_PF4.txt --mode cis -o ./Output_Sum_Final_chr_pos/Sub_fib
