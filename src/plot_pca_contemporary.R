library(tidyverse)
library(data.table)
library(viridis)

# read in metadata
sample_table <- fread("data/sample_table.csv")
sample_info <- sample_table[,c(1,2:4,6, 13:14)]

###################
#### LD PRUNED ####
###################

# read in data
pca <- fread("output/03_plink/ld_pruned_NZ_final//filtered_snps_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("output/03_plink/ld_pruned_NZ_final/filtered_snps_plink_pca.eigenval")

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

# longitude
ggplot(pca_info, aes(PC1, PC2, colour=NZ_longtitude_S))+
  geom_point(size = 3)+
  scale_colour_viridis()+ 
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# latitude
ggplot(pca_info, aes(PC1, PC2, colour=NZ_latitude_E))+
  geom_point(size = 3)+
  scale_colour_viridis()+ 
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


#######################
#### NOT LD PRUNED ####
#######################

# read in data
pca <- fread("output/03_plink/no_ldpruning_NZ_final//filtered_snps_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("output/03_plink/no_ldpruning_NZ_final/filtered_snps_plink_pca.eigenval")

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
long_pca_info <- melt(pca_info)

# plot all pcas - location
ggplot(long_pca_info, aes(sample_name, value, colour=Location))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)

# plot all pcas - east/west
ggplot(long_pca_info, aes(sample_name, value, colour=historic_E_W))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)