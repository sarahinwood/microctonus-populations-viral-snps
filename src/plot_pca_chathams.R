library(tidyverse)
library(data.table)
library(viridis)

# read in metadata
sample_table <- fread("data/sample_table.csv")
sample_info <- sample_table[,c(1:3,5:8,13,14,22)]

###################
#### LD PRUNED ####
###################

# read in data
pca <- fread("chathams_output/05_plink/ld_pruned_Chatham_PCA_75cov_only_filtered_SNPs/no_mind/filtered_SNPs_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("chathams_output/05_plink/ld_pruned_Chatham_PCA_75cov_only_filtered_SNPs/no_mind/filtered_SNPs_plink_pca.eigenval")

## what samples pruned out
sample_info_contemp_nz <- subset(sample_info, sample_timing=="contemporary")
pruned_out <- subset(sample_info_contemp_nz, !(sample_name %in% pca$sample_name))

# sort out the pca data
# remove extra sample name column
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample_name"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# merge with sample info
pca_info <- merge(pca, sample_info)
#pca_info$NZ_longtitude_S <- factor(pca_info$NZ_longtitude_S, levels=unique(pca_info$NZ_longtitude_S[order(pca_info$NZ_longtitude_S, decreasing = T)]))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+
  ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca - PC1 = LOCATION
ggplot(pca_info, aes(PC1, PC2, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# test other PC combos - only PC1 creates two subgroups of samples
ggplot(pca_info, aes(PC2, PC3, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))+
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

pca_info_to_melt <- pca_info[,-c(28, 29)]
long_pca_info <- melt(pca_info_to_melt)

long_pca_info$sample_name <- factor(long_pca_info$sample_name, levels=c("Murison_CI_4_pupa", "CI_x_B.H_26.8.21_1_pupa", "Chatham_Beach_Tuanui_1_pupa", "Chatham_Beach_Tuanui_2_pupa",   "Chatham_Beach_Tuanui_3_pupa", "Fortrose_1", "Fortrose_11", "Fortrose_14", "Fortrose_22", "Fortrose_23", "Fortrose_26", "Fortrose_4", "Lincoln_12", "Lincoln_13", "Lincoln_14", "Lincoln_4", "Lincoln_5", "Lincoln_Chatham_plate_18_live", "Lincoln_Chatham_plate_1_live", "Lincoln_Chatham_plate_2_live", "Lincoln_Chatham_plate_3_live", "Lincoln_Chatham_plate_4_live", "Lincoln_Chatham_plate_5_live", "Lincoln_Chatham_plate_9_live"))

# plot all pcas
ggplot(long_pca_info, aes(sample_name, value, colour=Location))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)

##what is DAPC - discriminant analysis of principle components?

#######################
#### NOT LD PRUNED ####
#######################

# read in data
pca <- fread("chathams_output/05_plink/no_ldpruning_Chatham_PCA_75cov_only_filtered_SNPs/no_mind/filtered_SNPs_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("chathams_output/05_plink/no_ldpruning_Chatham_PCA_75cov_only_filtered_SNPs/no_mind/filtered_SNPs_plink_pca.eigenval")

# what samples pruned out?
pruned_out_samples <- setdiff(sample_table$sample_name, pca$V1)
pruned_out_sample_table <- subset(sample_table, sample_name %in% pruned_out_samples)
pruned_out_sample_table <- subset(pruned_out_sample_table, sample_timing=="contemporary")

# sort out the pca data
# remove extra sample name column
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample_name"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# merge with sample info
pca_info <- merge(pca, sample_info)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+
  ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
ggplot(pca_info, aes(PC1, PC2, colour=Location))+ 
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# test other PC combos
ggplot(pca_info, aes(PC2, PC3, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))+
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

pca_info$sample_name <- factor(pca_info$sample_name, levels=pca_info$sample_name[order(pca_info$historic_E_W, pca_info$sample_name, decreasing=T)])
pca_info_to_melt <- pca_info[,-c(28, 29)]
long_pca_info <- melt(pca_info_to_melt)

# plot all pcas - location
ggplot(long_pca_info, aes(sample_name, value, colour=Location))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)

### ASW PCA

# read in data
pca <- fread("chathams_output/ASW_plink_no_ldpruning/mind_0.2/filtered_snps_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("chathams_output/ASW_plink_no_ldpruning/mind_0.2/filtered_snps_plink_pca.eigenval")

# sort out the pca data
# remove extra sample name column
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample_name"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

asw_sample_info <- fread("../../asw_projects/asw-pangenome/data/sample_table.csv")

# merge with sample info
pca_info <- merge(pca, asw_sample_info)

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+
  ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
ggplot(pca_info, aes(PC1, PC2, colour=Location))+ 
  geom_point(size = 3, alpha=0.7)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
