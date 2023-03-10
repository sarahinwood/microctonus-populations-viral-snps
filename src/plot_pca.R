library(tidyverse)
library(data.table)
library(viridis)

# read in metadata
sample_table <- fread("data/sample_table.csv")
sample_info <- sample_table[,c(1,6,8,9)]

#######################
#### NOT LD PRUNED ####
#######################

# read in data
pca <- fread("output/03_plink/no_ldpruning_final/filtered_snps_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("output/03_plink/no_ldpruning_final/filtered_snps_plink_pca.eigenval")

# what samples pruned out?
pruned_out_samples <- setdiff(sample_table$sample_name, pca$V1)
pruned_out_sample_table <- subset(sample_table, sample_name %in% pruned_out_samples)

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
ggplot(pca_info, aes(PC1, PC2, colour=Location, shape=sample_timing))+ 
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  labs(shape="Sample timing")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# test other PC combos
ggplot(pca_info, aes(PC2, PC3, shape=sample_timing, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  labs(shape="Sample timing")+
  xlab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))+
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

pca_info$sample_name <- factor(pca_info$sample_name, levels=pca_info$sample_name[order(pca_info$sample_timing, pca_info$sample_name, decreasing=T)])
long_pca_info <- melt(pca_info)

# plot all pcas - location
ggplot(long_pca_info, aes(sample_name, value, colour=Location, shape=sample_timing))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  labs(shape="Sample timing")+
  facet_wrap(~variable)

# plot all pcas - east/west
ggplot(long_pca_info, aes(sample_name, value, colour=sample_timing))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)

### only historic ###
pca_info_historic <- subset(pca_info, sample_timing=="historic")
ggplot(pca_info_historic, aes(PC1, PC2, shape=historic_E_W, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
##########

###################
#### LD PRUNED ####
###################

# read in data
pca <- fread("output/03_plink/ld_pruned/filtered_snps_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("output/03_plink/ld_pruned/filtered_snps_plink_pca.eigenval")

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

# plot pca - PC1 = LOCATION
ggplot(pca_info, aes(PC1, PC2, shape=sample_timing, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

### only historic ###
pca_info_historic <- subset(pca_info, sample_timing=="historic")
ggplot(pca_info_historic, aes(PC1, PC2, shape=historic_E_W, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# test other PC combos - only PC1 creates two subgroups of samples
ggplot(pca_info_historic, aes(PC2, PC3, shape=historic_E_W, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))+
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

long_pca_info <- melt(pca_info)

# plot all pcas
ggplot(long_pca_info, aes(sample_name, value, colour=Location))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)

##what is DAPC - discriminant analysis of principle components?
