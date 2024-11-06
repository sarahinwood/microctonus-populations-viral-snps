library(tidyverse)
library(data.table)
library(viridis)

# read in metadata
sample_table <- fread("data/sample_table.csv")
sample_table_histSA <- subset(sample_table, sample_timing=="historic rearing")
sample_info <- sample_table_histSA[,c(1:3,5:7,13,14,22)]

###################
#### LD PRUNED ####
###################

# read in data
pca <- fread("output/05_plink/ld_pruned_SouthAmerican_Mhyp_75cov_only/mind_0.2/filtered_SNPs_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("output/05_plink/ld_pruned_SouthAmerican_Mhyp_75cov_only/mind_0.2/filtered_SNPs_plink_pca.eigenval")

# what samples pruned out?
pruned_out_samples <- setdiff(sample_table_histSA$sample_name, pca$V1)
pruned_out_sample_table <- subset(sample_table, sample_name %in% pruned_out_samples)

# sort out the pca data
# remove extra sample name column
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample_name"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# merge with sample info
pca_info <- merge(pca, sample_info, by="sample_name")
pca_info$Location <- factor(pca_info$Location, levels=c("Concepcion (South Chile)", "La Serena",
                                                        "San Carlos de Bariloche",
                                                        "Mendoza", "Rio Negro", "Hilario Ascasubi", "Colonia (Uruguay)", "Porto Alegre"))
pca_info$historic_E_W <- factor(pca_info$historic_E_W, levels=c("West", "Intermediate", "East"))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+
  ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
ggplot(pca_info, aes(PC1, PC2, colour=Location, shape=historic_E_W))+ 
  geom_point(size = 3, alpha=0.8)+
  scale_colour_manual(values=c("Concepcion (South Chile)" = "#8c2981",
                               "La Serena" = "#de4968",
                               "San Carlos de Bariloche" = "#fe9f6d",
                               "Mendoza" = "#414487",
                               "Rio Negro" = "#2a788e",
                               "Hilario Ascasubi" = "#22a884",
                               "Colonia (Uruguay)" = "#7ad151",
                               "Porto Alegre" = "#fde725"))+
  scale_shape_manual(values=c("West"=15,
                              "Intermediate"=17,
                              "East" = 19))+
  coord_equal()+
  theme_light()+
  labs(shape="East vs West of Andes")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

pca_info$rearing_timing <- ifelse(pca_info$seq_cohort=="cohort4", paste("Early rearing"), paste("Late rearing"))
# plot pca
ggplot(pca_info, aes(PC1, PC2, colour=Location, shape=rearing_timing))+ 
  geom_point(size = 3, alpha=0.8)+
  scale_colour_manual(values=c("Concepcion (South Chile)" = "#8c2981",
                               "La Serena" = "#de4968",
                               "San Carlos de Bariloche" = "#fe9f6d",
                               "Mendoza" = "#414487",
                               "Rio Negro" = "#2a788e",
                               "Hilario Ascasubi" = "#22a884",
                               "Colonia (Uruguay)" = "#7ad151",
                               "Porto Alegre" = "#fde725"))+
  coord_equal()+
  theme_light()+
  labs(shape="Rearing stage")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# coloured by seq batch
pca_info$sample_batch <- ifelse(pca_info$seq_cohort=="cohort4", paste("Batch 2 (earlier)"), paste("Batch 1 (later)"))
ggplot(pca_info, aes(PC1, PC2, colour=sample_batch, shape=historic_E_W))+ 
  geom_point(size = 3, alpha=0.8)+
  scale_colour_viridis(discrete=T)+
  scale_shape_manual(values=c("West"=15,
                              "Intermediate"=17,
                              "East" = 19))+
  coord_equal()+
  theme_light()+
  labs(shape="East vs West of Andes", colour="Sample Batch")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))


ggplot(pca_info, aes(PC1, PC2, colour=historic_E_W))+ 
  geom_point(size = 3, alpha=0.9)+
  scale_colour_viridis(discrete=T)+
  coord_equal()+
  theme_light()+
  labs(shape="East vs West of Andes")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## plot pca longitude
ggplot(pca_info, aes(PC1, PC2, colour=historic_longitude))+
  geom_point(size = 3, alpha=0.9)+
  scale_colour_viridis()+
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## plot pca latitude
ggplot(pca_info, aes(PC1, PC2, colour=historic_latitude))+
  geom_point(size = 3, alpha=0.9)+
  scale_colour_viridis()+
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# test other PC combos
ggplot(pca_info, aes(PC2, PC3, colour=historic_E_W))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  labs(shape="East vs West of Andes")+
  xlab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))+
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

pca_info$sample_name <- factor(pca_info$sample_name, levels=pca_info$sample_name[order(pca_info$historic_E_W, pca_info$sample_name, decreasing=T)])
pca_info_melt <- pca_info[,-c(24, 26, 27)]
long_pca_info <- melt(pca_info_melt)
long_pca_info$historic_longitude <- as.numeric(long_pca_info$historic_longitude)

# plot all pcas - location
ggplot(long_pca_info, aes(sample_name, value, colour=Location, shape=historic_E_W))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  labs(shape="East vs West of Andes")+
  facet_wrap(~variable)

# plot all pcas - east/west
ggplot(long_pca_info, aes(sample_name, value, colour=historic_E_W))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)


pca_info_lat <- pca_info
pca_info_lat$historic_longitude <- as.numeric(pca_info_lat$historic_longitude)
pca_info_lat$sample_name <- factor(pca_info_lat$sample_name, levels=pca_info_lat$sample_name[order(pca_info_lat$historic_longitude, pca_info_lat$sample_name, decreasing=T)])
sample_to_lat <- pca_info_lat[,c(1,25)]
pca_info_lat <- pca_info_lat[,-25]
long_pca_info_lat <- melt(pca_info_lat)
long_pca_info_lat <- merge(long_pca_info_lat, sample_to_lat)

# plot all pcas - longitude
ggplot(long_pca_info_lat, aes(sample_name, value, colour=historic_longitude))+
  geom_point()+
  scale_colour_viridis()+ 
  theme_light()+
  facet_wrap(~variable)


#######################
#### NOT LD PRUNED ####
#######################

# read in data
pca <- fread("output/05_plink/no_ldpruning_SouthAmerican_Mhyp_75cov_only/mind_0.2/filtered_SNPs_plink_pca.eigenvec", header = FALSE)
eigenval <- scan("output/05_plink/no_ldpruning_SouthAmerican_Mhyp_75cov_only/mind_0.2/filtered_SNPs_plink_pca.eigenval")

# what samples pruned out?
pruned_out_samples <- setdiff(sample_table_histSA$sample_name, pca$V1)
pruned_out_sample_table <- subset(sample_table_histSA, sample_name %in% pruned_out_samples)

# sort out the pca data
# remove extra sample name column
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample_name"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# merge with sample info
pca_info <- merge(pca, sample_info, by="sample_name")
pca_info$Location <- factor(pca_info$Location, levels=c("Concepcion (South Chile)", "La Serena",
                                                        "San Carlos de Bariloche",
                                                        "Mendoza", "Rio Negro", "Hilario Ascasubi", "Colonia (Uruguay)", "Porto Alegre"))
pca_info$historic_E_W <- factor(pca_info$historic_E_W, levels=c("West", "Intermediate", "East"))


# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")+
  ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

pca_info$rearing_timing <- ifelse(pca_info$seq_cohort=="cohort4", paste("Early rearing"), paste("Late rearing"))
# plot pca
ggplot(pca_info, aes(PC1, PC2, colour=Location, shape=rearing_timing))+ 
  geom_point(size = 3, alpha=0.8)+
  scale_colour_manual(values=c("Concepcion (South Chile)" = "#8c2981",
                               "La Serena" = "#de4968",
                               "San Carlos de Bariloche" = "#fe9f6d",
                               "Mendoza" = "#414487",
                               "Rio Negro" = "#2a788e",
                               "Hilario Ascasubi" = "#22a884",
                               "Colonia (Uruguay)" = "#7ad151",
                               "Porto Alegre" = "#fde725"))+
  coord_equal()+
  theme_light()+
  labs(shape="Rearing stage")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## same plot only early
pca_info_early <- subset(pca_info, rearing_timing=="Early rearing")
ggplot(pca_info_early, aes(PC1, PC2, colour=Location, shape=rearing_timing))+ 
  geom_point(size = 3, alpha=0.8)+
  scale_colour_manual(values=c("Concepcion (South Chile)" = "#8c2981",
                               "La Serena" = "#de4968",
                               "San Carlos de Bariloche" = "#fe9f6d",
                               "Mendoza" = "#414487",
                               "Rio Negro" = "#2a788e",
                               "Hilario Ascasubi" = "#22a884",
                               "Colonia (Uruguay)" = "#7ad151",
                               "Porto Alegre" = "#fde725"))+
  coord_equal()+
  theme_light()+
  labs(shape="Rearing stage")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))


# plot pca
ggplot(pca_info, aes(PC1, PC2, colour=Location, shape=historic_E_W))+ 
  geom_point(size = 3, alpha=0.9)+
  scale_colour_viridis(discrete=T, option="D")+
  coord_equal()+
  theme_light()+
  labs(shape="East vs West of Andes")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## plot pca longitude
ggplot(pca_info, aes(PC1, PC2, colour=historic_longitude))+
  geom_point(size = 3, alpha=0.9)+
  scale_colour_viridis()+
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

## plot pca longitude
ggplot(pca_info, aes(PC1, PC2, colour=historic_latitude))+
  geom_point(size = 3, alpha=0.9)+
  scale_colour_viridis()+
  coord_equal()+
  theme_light()+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# test other PC combos
ggplot(pca_info, aes(PC2, PC3, shape=historic_E_W, colour=Location))+
  geom_point(size = 3)+
  scale_colour_viridis(discrete=T)+ 
  coord_equal()+
  theme_light()+
  labs(shape="East vs West of Andes")+
  xlab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))+
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

pca_info$sample_name <- factor(pca_info$sample_name, levels=pca_info$sample_name[order(pca_info$historic_E_W, pca_info$sample_name, decreasing=T)])
pca_info$historic_longitude <- as.character(pca_info$historic_longitude)
pca_info_melt <- pca_info[,-c(24, 26, 27)]
long_pca_info <- melt(pca_info_melt)
long_pca_info$historic_longitude <- as.numeric(long_pca_info$historic_longitude)

# plot all pcas - location
ggplot(long_pca_info, aes(sample_name, value, colour=Location, shape=historic_E_W))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  labs(shape="East vs West of Andes")+
  facet_wrap(~variable)

# plot all pcas - east/west
ggplot(long_pca_info, aes(sample_name, value, colour=historic_E_W))+
  geom_point()+
  scale_colour_viridis(discrete=T)+ 
  theme_light()+
  facet_wrap(~variable)

# plot all pcas - longitude
ggplot(long_pca_info, aes(sample_name, value, colour=historic_longitude))+
  geom_point()+
  scale_colour_viridis()+ 
  theme_light()+
  facet_wrap(~variable)