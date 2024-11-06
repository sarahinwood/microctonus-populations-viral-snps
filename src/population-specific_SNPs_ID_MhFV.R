library(data.table)
library(tidyverse)

##############
## Ascasubi ##
##############

Ascasubi_allele_freqs <- data.table(read_tsv("output/05_plink/population_specific_SNPs/ecotype/ecotype.Ascasubi.afreq", col_names=T, id="path"))
Ascasubi_allele_freqs_dt <- Ascasubi_allele_freqs[,c(3:6)]
setnames(Ascasubi_allele_freqs_dt, old="ALT_FREQS", new="Ascasubi_ALT_FREQS")

##############
## Mendoza ##
##############

Mendoza_allele_freqs <- data.table(read_tsv("output/05_plink/population_specific_SNPs/ecotype/ecotype.Mendoza.afreq", col_names=T, id="path"))
Mendoza_allele_freqs_dt <- Mendoza_allele_freqs[,c(3:6)]
setnames(Mendoza_allele_freqs_dt, old="ALT_FREQS", new="Mendoza_ALT_FREQS")

##############
## RioNegro ##
##############

RioNegro_allele_freqs <- data.table(read_tsv("output/05_plink/population_specific_SNPs/ecotype/ecotype.RioNegro.afreq", col_names=T, id="path"))
RioNegro_allele_freqs_dt <- RioNegro_allele_freqs[,c(3:6)]
setnames(RioNegro_allele_freqs_dt, old="ALT_FREQS", new="RioNegro_ALT_FREQS")

###########
## MERGE ##
###########

ecotype_freqs <- left_join(Ascasubi_allele_freqs_dt, left_join(Mendoza_allele_freqs_dt, RioNegro_allele_freqs_dt))
ecotype_freqs[ecotype_freqs == "NaN"] <- NA

specific_SNPs_Mendoza <- subset(ecotype_freqs, Mendoza_ALT_FREQS==1&Ascasubi_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Mendoza$ecotype_specificity <- paste("Mendoza")
#format to write specific SNP list
specific_SNPs_Mendoza$scaffold <- tstrsplit(specific_SNPs_Mendoza$ID, ':', keep=1)
specific_SNPs_Mendoza$bp <- tstrsplit(specific_SNPs_Mendoza$ID, ':', keep=2)
Mendoza_SNP_list <- specific_SNPs_Mendoza[,c(8,9)]
write_tsv(Mendoza_SNP_list, "output/05_plink/population_specific_SNPs/ecotype/Mendoza_specific_SNPs.tsv")

specific_SNPs_Ascasubi <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==1&Mendoza_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Ascasubi$ecotype_specificity <- paste("Ascasubi")
#format to write specific SNP list
specific_SNPs_Ascasubi$scaffold <- tstrsplit(specific_SNPs_Ascasubi$ID, ':', keep=1)
specific_SNPs_Ascasubi$bp <- tstrsplit(specific_SNPs_Ascasubi$ID, ':', keep=2)
Ascasubi_SNP_list <- specific_SNPs_Ascasubi[,c(8,9)]
write_tsv(Ascasubi_SNP_list, "output/05_plink/population_specific_SNPs/ecotype/Ascasubi_specific_SNPs.txt", col_names=F)

specific_SNPs_RioNegro <- subset(ecotype_freqs, RioNegro_ALT_FREQS==1&Ascasubi_ALT_FREQS==0&Mendoza_ALT_FREQS==0)
specific_SNPs_RioNegro$ecotype_specificity <- paste("RioNegro")
#format to write specific SNP list
specific_SNPs_RioNegro$scaffold <- tstrsplit(specific_SNPs_RioNegro$ID, ':', keep=1)
specific_SNPs_RioNegro$bp <- tstrsplit(specific_SNPs_RioNegro$ID, ':', keep=2)
RioNegro_SNP_list <- specific_SNPs_RioNegro[,c(8,9)]
write_tsv(RioNegro_SNP_list, "output/05_plink/population_specific_SNPs/ecotype/RioNegro_specific_SNPs.tsv")


ecotype_specific_SNPs <- full_join(specific_SNPs_Mendoza, full_join(specific_SNPs_Ascasubi, specific_SNPs_RioNegro))
ecotype_specific_SNP_list <- ecotype_specific_SNPs[,c(1,7)]

#####################################
## NZ samples - early post-release ##
#####################################

HawkesBay_allele_freqs <- data.table(read_tsv("output/05_plink/population_specific_SNPs/NZ_postrelease_early/location.HawkesBay.afreq", col_names=T, id="path"))
HawkesBay_allele_freqs_dt <- HawkesBay_allele_freqs[,c(3:6)]
setnames(HawkesBay_allele_freqs_dt, old="ALT_FREQS", new="HawkesBay_ALT_FREQS")

Reporoa_allele_freqs <- data.table(read_tsv("output/05_plink/population_specific_SNPs/NZ_postrelease_early/location.Reporoa.afreq", col_names=T, id="path"))
Reporoa_allele_freqs_dt <- Reporoa_allele_freqs[,c(3:6)]
setnames(Reporoa_allele_freqs_dt, old="ALT_FREQS", new="Reporoa_ALT_FREQS")

Ruakura_allele_freqs <- data.table(read_tsv("output/05_plink/population_specific_SNPs/NZ_postrelease_early/location.Ruakura.afreq", col_names=T, id="path"))
Ruakura_allele_freqs_dt <- Ruakura_allele_freqs[,c(3:6)]
setnames(Ruakura_allele_freqs_dt, old="ALT_FREQS", new="Ruakura_ALT_FREQS")

Wellsford_allele_freqs <- data.table(read_tsv("output/05_plink/population_specific_SNPs/NZ_postrelease_early/location.Wellsford.afreq", col_names=T, id="path"))
Wellsford_allele_freqs_dt <- Wellsford_allele_freqs[,c(3:6)]
setnames(Wellsford_allele_freqs_dt, old="ALT_FREQS", new="Wellsford_ALT_FREQS")

NZ_early_allele_freqs <- left_join(HawkesBay_allele_freqs_dt, left_join(Reporoa_allele_freqs_dt, left_join(Ruakura_allele_freqs_dt, Wellsford_allele_freqs_dt)))
NZ_early_allele_freqs[NZ_early_allele_freqs == "NaN"] <- NA

NZ_early_freqs_ecotype_specificity <- merge(ecotype_specific_SNP_list, NZ_early_allele_freqs)
