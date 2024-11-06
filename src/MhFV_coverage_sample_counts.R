library(tidyverse)
library(data.table)

sample_table <- fread("data/sample_table.csv")

sample_table_85 <- subset(sample_table, MhFV_coverage>85)
sample_table_75 <- subset(sample_table, MhFV_coverage>75)
sample_table_70 <- subset(sample_table, MhFV_coverage>70)

sample_table_85_counts <- count(sample_table_85, Location, sample_timing)
sample_table_75_counts <- count(sample_table_75, Location, sample_timing)
sample_table_70_counts <- count(sample_table_70, Location, sample_timing)

setnames(sample_table_85_counts, old="n", new="n_85cov")
setnames(sample_table_75_counts, old="n", new="n_75cov")
setnames(sample_table_70_counts, old="n", new="n_70cov")

cov_sample_counts <- merge(sample_table_85_counts, merge(sample_table_75_counts, sample_table_70_counts))

cov_sample_counts$diff_75_85 <- (cov_sample_counts$n_75cov-cov_sample_counts$n_85cov)
cov_sample_counts$diff_70_75 <- (cov_sample_counts$n_70cov-cov_sample_counts$n_75cov)





ggplot(sample_table_histSA, aes(Location, MhFV_coverage, colour=label))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2, alpha = 0.75,shape = 16)+
  theme_bw()+
  xlab("Sample Location")+ylab("MhFV genome sequencing coverage")+labs(colour="Sample type")
