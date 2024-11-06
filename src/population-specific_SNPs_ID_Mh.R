library(data.table)
library(tidyverse)

## how does diploidy affect this approach - when freq = 1 does that mean present or homozygous?

## all heterozygous

##############
## Ascasubi ##
##############

Ascasubi_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.Ascasubi.afreq", col_names=T, id="path"))
Ascasubi_allele_freqs <- Ascasubi_allele_freqs[,c(3:6)]
setnames(Ascasubi_allele_freqs, old="ALT_FREQS", new="Ascasubi_ALT_FREQS")

###############
## Bariloche ##
###############

Bariloche_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.Bariloche.afreq", col_names=T, id="path"))
Bariloche_allele_freqs <- Bariloche_allele_freqs[,c(3:6)]
setnames(Bariloche_allele_freqs, old="ALT_FREQS", new="Bariloche_ALT_FREQS")

#############
## Colonia ##
#############

Colonia_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.Colonia.afreq", col_names=T, id="path"))
Colonia_allele_freqs <- Colonia_allele_freqs[,c(3:6)]
setnames(Colonia_allele_freqs, old="ALT_FREQS", new="Colonia_ALT_FREQS")

################
## Concepcion ##
################

Concepcion_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.Concepcion.afreq", col_names=T, id="path"))
Concepcion_allele_freqs <- Concepcion_allele_freqs[,c(3:6)]
setnames(Concepcion_allele_freqs, old="ALT_FREQS", new="Concepcion_ALT_FREQS")

##############
## LaSerena ##
##############

LaSerena_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.LaSerena.afreq", col_names=T, id="path"))
LaSerena_allele_freqs <- LaSerena_allele_freqs[,c(3:6)]
setnames(LaSerena_allele_freqs, old="ALT_FREQS", new="LaSerena_ALT_FREQS")

#############
## Mendoza ##
#############

Mendoza_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.Mendoza.afreq", col_names=T, id="path"))
Mendoza_allele_freqs <- Mendoza_allele_freqs[,c(3:6)]
setnames(Mendoza_allele_freqs, old="ALT_FREQS", new="Mendoza_ALT_FREQS")

#################
## PortoAlegre ##
#################

PortoAlegre_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.PortoAlegre.afreq", col_names=T, id="path"))
PortoAlegre_allele_freqs <- PortoAlegre_allele_freqs[,c(3:6)]
setnames(PortoAlegre_allele_freqs, old="ALT_FREQS", new="PortoAlegre_ALT_FREQS")

##############
## RioNegro ##
##############

RioNegro_allele_freqs <- data.table(read_tsv("mh_private_alleles/ecotype_specific_SNPs/ecotype.RioNegro.afreq", col_names=T, id="path"))
RioNegro_allele_freqs <- RioNegro_allele_freqs[,c(3:6)]
setnames(RioNegro_allele_freqs, old="ALT_FREQS", new="RioNegro_ALT_FREQS")

###########
## MERGE ## set up with a min freq above 0 e.g. 0.3
###########

ecotype_freqs <- left_join(Ascasubi_allele_freqs,
                           left_join(Bariloche_allele_freqs,
                                     left_join(Colonia_allele_freqs,
                                               left_join(Concepcion_allele_freqs,
                                                         left_join(LaSerena_allele_freqs,
                                                                   left_join(Mendoza_allele_freqs,
                                                                             left_join(PortoAlegre_allele_freqs, RioNegro_allele_freqs)))))))
ecotype_freqs[ecotype_freqs == "NaN"] <- NA

specific_SNPs_Ascasubi <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==1&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Ascasubi$ecotype_specificity <- paste("Ascasubi")

specific_SNPs_Ascasubi_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS>0.1&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Ascasubi_not1$ecotype_specificity <- paste("Ascasubi")

specific_SNPs_Bariloche <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==1&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Bariloche$ecotype_specificity <- paste("Bariloche")

specific_SNPs_Bariloche_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS>0.1&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Bariloche_not1$ecotype_specificity <- paste("Bariloche")

specific_SNPs_Colonia <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==1&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Colonia$ecotype_specificity <- paste("Colonia")

specific_SNPs_Colonia_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS>0.1&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Colonia_not1$ecotype_specificity <- paste("Colonia")

specific_SNPs_Concepcion <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==1&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Concepcion$ecotype_specificity <- paste("Concepcion")

specific_SNPs_Concepcion_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS>0.1&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Concepcion_not1$ecotype_specificity <- paste("Concepcion")

specific_SNPs_LaSerena <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==1&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==00)
specific_SNPs_LaSerena$ecotype_specificity <- paste("LaSerena")

specific_SNPs_LaSerena_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS>0.1&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==00)
specific_SNPs_LaSerena_not1$ecotype_specificity <- paste("LaSerena")

specific_SNPs_Mendoza <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==1&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Mendoza$ecotype_specificity <- paste("Mendoza")

specific_SNPs_Mendoza_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS>0.1&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_Mendoza_not1$ecotype_specificity <- paste("Mendoza")

specific_SNPs_PortoAlegre <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==1&RioNegro_ALT_FREQS==0)
specific_SNPs_PortoAlegre$ecotype_specificity <- paste("PortoAlegre")

specific_SNPs_PortoAlegre_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS>0.1&RioNegro_ALT_FREQS==0)
specific_SNPs_PortoAlegre_not1$ecotype_specificity <- paste("PortoAlegre")

specific_SNPs_RioNegro <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==1)
specific_SNPs_RioNegro$ecotype_specificity <- paste("RioNegro")

specific_SNPs_RioNegro_not1 <- subset(ecotype_freqs, Ascasubi_ALT_FREQS==0&Bariloche_ALT_FREQS==0&Colonia_ALT_FREQS==0&Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS>0.1)
specific_SNPs_RioNegro_not1$ecotype_specificity <- paste("RioNegro")

ecotype_specific_SNPs_not1 <- full_join(specific_SNPs_Mendoza_not1,
                                        full_join(specific_SNPs_Ascasubi_not1,
                                                  full_join(specific_SNPs_Bariloche_not1,
                                                            full_join(specific_SNPs_Colonia_not1,
                                                                      full_join(specific_SNPs_Concepcion_not1,
                                                                                full_join(specific_SNPs_LaSerena_not1,
                                                                                          full_join(specific_SNPs_Mendoza_not1,
                                                                                                    full_join(specific_SNPs_PortoAlegre_not1, specific_SNPs_RioNegro_not1))))))))
ecotype_specific_SNPs_list_not1 <- ecotype_specific_SNPs_not1[,c(1,12)]
ecotype_specific_SNPs_list_not1$scaffold <- tstrsplit(ecotype_specific_SNPs_list_not1$ID, "_", keep=2)
ecotype_specific_SNPs_list_not1$scaffold <- paste("scaffold", ecotype_specific_SNPs_list_not1$scaffold, sep="_")
ecotype_specific_SNPs_list_not1$bp <- tstrsplit(ecotype_specific_SNPs_list_not1$ID, "_", keep=3)
ecotype_specific_SNPs_list_sites <- ecotype_specific_SNPs_list_not1[,c(3,4)]
## filter VCF to ecotype SNPs
write_tsv(ecotype_specific_SNPs_list_sites, "mh_private_alleles/ecotype_SNPs_above_freq_0.1.txt", col_names=F)



ecotype_counts <- ecotype_specific_SNPs_list_not1 %>% group_by(ecotype_specificity) %>% summarise(count = n())
western <- c("LaSerena", "Concepcion")
intermediate <- c("Bariloche")
ecotype_counts$east_west <- ifelse(ecotype_counts$ecotype_specificity=="LaSerena", paste("western"),
                                   ifelse(ecotype_counts$ecotype_specificity=="Concepcion", paste("western"),
                                          ifelse(ecotype_counts$ecotype_specificity=="Bariloche", paste("intermediate"), paste("eastern"))))


## I think plink is counting SNPs missing calls as not being REF so adding to ALT counts?
## do for mitochondrial also

##################
## East vs West ##
##################

# east or west specific and freq = 1

specific_SNPs_West <- subset(ecotype_freqs, Concepcion_ALT_FREQS==1&LaSerena_ALT_FREQS==1&Ascasubi_ALT_FREQS==0&Colonia_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_West$ecotype_specificity <- paste("West")

specific_SNPs_East <- subset(ecotype_freqs, Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Ascasubi_ALT_FREQS==1&Colonia_ALT_FREQS==1&Mendoza_ALT_FREQS==1&PortoAlegre_ALT_FREQS==1&RioNegro_ALT_FREQS==1)
specific_SNPs_East$ecotype_specificity <- paste("East")

# east or west specific but freq is not 1

specific_SNPs_West_not1 <- subset(ecotype_freqs, Concepcion_ALT_FREQS!=0&LaSerena_ALT_FREQS!=0&Ascasubi_ALT_FREQS==0&Colonia_ALT_FREQS==0&Mendoza_ALT_FREQS==0&PortoAlegre_ALT_FREQS==0&RioNegro_ALT_FREQS==0)
specific_SNPs_West_not1$ecotype_specificity <- paste("West")

specific_SNPs_East_not1 <- subset(ecotype_freqs, Concepcion_ALT_FREQS==0&LaSerena_ALT_FREQS==0&Ascasubi_ALT_FREQS!=0&Colonia_ALT_FREQS!=0&Mendoza_ALT_FREQS!=0&PortoAlegre_ALT_FREQS!=0&RioNegro_ALT_FREQS!=0)
specific_SNPs_East_not1$ecotype_specificity <- paste("East")


west_east_specific_SNPs_not1 <- full_join(specific_SNPs_West_not1, specific_SNPs_East_not1)
west_east_specific_SNPs_list_not1 <- west_east_specific_SNPs_not1[,c(1,12)]

## filter VCF to east/west SNPs
specific_SNP_list <- full_join(west_east_specific_SNPs_list_not1, ecotype_specific_SNPs_list_not1)
fwrite(list(west_east_specific_SNPs_list_not1$ID), "mh_private_alleles/east_west_specific_SNPs.txt")

