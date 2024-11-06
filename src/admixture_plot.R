library(tidyverse)
library(ggplot2)
library(forcats)
library(data.table)
library(viridis)

## which K value best?
cv_res <- fread('output/05_plink/ld_pruned_SouthAmericanCohort4_75cov_only/mind_0.2/admixture_cv_res.out', header=F)
setorder(cv_res, V4)

## read in admixture for best K value
admixture_9 <- read_delim("output/100_plink/admixture_10_10_0.2/admixture.9.Q",
           col_names = paste0("Q",seq(1:3)),
           delim=" ")
admixture_9$k <- "9"

sample_names <- fread("output/05_plink/ld_pruned_SouthAmericanCohort4_Mhyp_75cov_only/mind_0.2/filtered_SNPs_plink_pca", header=F)
admixture_9$sample <- sample_names$V1


long_admix3 <- gather(admixture_3, Q, value, -sample,-k)
## get location names
long_admix3 <- data.table(long_admix3)
long_admix3$Location <- tstrsplit(long_admix3$sample, "geo_", keep=2)
long_admix3$Location <- gsub("[[:digit:]]","", long_admix3$Location)
long_admix3$Location <- gsub("O-", "Ophir", long_admix3$Location)
long_admix3$Location <- gsub("-", " ", long_admix3$Location)
long_admix3$Location <- ifelse(grepl("Chatham", long_admix3$sample), paste("Chatham Island"),
                                     ifelse(grepl("Stewart", long_admix3$sample), paste("Stewart Island"),
                                                  ifelse(grepl("Lincoln", long_admix3$sample), paste("Lincoln"),
                                                               ifelse(grepl("Fortrose", long_admix3$sample), paste("Fortrose"), paste(long_admix3$Location)))))

north_island <- c("Coromandel", "Ruakura", "Taranaki", "Wellington")
long_admix3$Island <- ifelse(long_admix3$Location %in% north_island, paste("North Island"),
                          ifelse(long_admix3$Location=="Chatham Island", paste("Chatham Island"),
                                 ifelse(long_admix3$Location=="Stewart Island", paste("Stewart Island"), paste("South Island"))))

north_main_divide <- c("Coromandel", "Ruakura", "Taranaki", "Wellington", "Greymouth")
long_admix3$main_divide <- ifelse(long_admix3$Location %in% north_main_divide, paste("North"),
                               ifelse(long_admix3$Location=="Chatham Island", paste("Chatham Island"),
                                      ifelse(long_admix3$Location=="Stewart Island", paste("Stewart Island"),
                                             paste("South"))))

long_admix3$main_divide <- factor(long_admix3$main_divide, levels=c("North", "South", "Stewart Island", "Chatham Island"))
long_admix3$Location <- factor(long_admix3$Location, levels=c("Coromandel", "Ruakura", "Taranaki", "Wellington", "Reefton", "Greymouth", "Lincoln", "Ophir", "Mararoa Downs", "Mossburn", "Fortrose", "Stewart Island", "Chatham Island"))


sample_order <- c()
long_admix3$sample <- factor(long_admix3$sample, levels=sample_order)

ggplot(long_admix3, aes(x=sample, y=value, fill=factor(Q)))+ 
  geom_bar(stat="identity",position="stack")+
  xlab("Sample")+
  ylab("Ancestry")+
  theme_bw()+
  scale_fill_viridis(discrete=T)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

