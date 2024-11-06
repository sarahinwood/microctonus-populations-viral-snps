library(tidyverse)
library(data.table)
library(viridis)

sample_table <- fread("data/sample_table.csv")
sample_table_histSA <- subset(sample_table, sample_timing=="historic rearing")
sample_info <- sample_table_histSA[,c(1:3,5:7,13,14,22)]


sample_info$Location <- factor(sample_info$Location, levels=c("Concepcion (South Chile)", "La Serena",
                                                                                    "San Carlos de Bariloche",
                                                                                    "Mendoza", "Rio Negro", "Hilario Ascasubi", "Colonia (Uruguay)", "Porto Alegre"))
sample_info$historic_E_W <- factor(sample_info$historic_E_W, levels=c("West", "Intermediate", "East"))
sample_info$rearing_timing <- ifelse(sample_info$seq_cohort=="cohort4", paste("Early rearing"), paste("Late rearing"))

sample_info$plot_group <- paste(sample_info$Location, sample_info$rearing_timing)

ggplot(sample_info, aes(Location, MhFV_coverage, color=rearing_timing))+
  geom_boxplot(fill = NA,
               aes(group = plot_group, colour=rearing_timing),
               outlier.size = -1, # prevents plotting of outliers which already plotted by geom_point
               show.legend = FALSE)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.75),
             alpha=0.7,
             size=2.5,
             shape=16)+
  scale_colour_manual(values=c("Early rearing"="#21918c",
                               "Late rearing"="#440154"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust=1))+
  labs(colour="Rearing stage")+
  ylab("Sequence coverage of MhFV genome (%)")+
  xlab("Ecotype collection location")



pca_info_c4 <- subset(pca_info, seq_cohort=="cohort4")
pca_info_c2 <- subset(pca_info, seq_cohort!="cohort4")


sample_table <- fread("data/sample_table.csv")
#sample_table_c4 <- subset(sample_table, seq_cohort=="cohort4")
sample_table_SA <- subset(sample_table, sample_timing=="historic rearing")
sample_table_SA$label <- ifelse(sample_table_SA$seq_cohort=="cohort4", paste("Eppendorf"), paste("Plate"))
sample_table_SA$xlabel <- paste(sample_table_SA$Location, sample_table_SA$label, sep=" ")

ggplot(sample_table_SA, aes(xlabel, MhFV_coverage, colour=label))+
  geom_boxplot()+
  geom_point(size = 2, alpha = 0.75,shape = 16)+
  theme_bw()+
  xlab("Sample Location")+ylab("MhFV genome sequencing coverage")+labs(colour="Sample type")
