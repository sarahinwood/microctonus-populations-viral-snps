library(data.table)
library(tidyverse)
library(ggplot2)

##https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

# FS  - is the Phred-scaled probability that there is strand bias at that site, 0 = little to no bias, GATK rec.s filtering at least 60
# QD - quality by depth, quality score normalised by read depth, GATK rec.s remove those 2.0 or less
# DP > 10 - depth above 10
# DP = sum of read depth over all samples, Variants that have a high mean depth across all samples are indicative of either repetitive regions or paralogs. Remove those with twice mean
# SOR - measurement for strand bias, but it takes into account the ratios of reads that cover both alleles
# MQ - mapping quality, good = ~60
# ReadPosRankSum - measures the relative position of the reference vs the alternative allele within reads. E.g., SNPs towards the ends of reads may be more likely a sequencing error, < -8.0

unfSnps_variant_stats <- fread("output/02_filtering/unfiltered_SNPs_stats.txt", header=F, na.strings=".")
unfSnps_variant_stats$filtration_status <- paste("unfiltered")

filSnps_variant_stats <- fread("output/02_filtering/pass-filtered_SNPs_stats.txt", header=F, na.strings=".")
filSnps_variant_stats$filtration_status <- paste("filtered")

#join tables
variant_stats <- full_join(unfSnps_variant_stats, filSnps_variant_stats)
setnames(variant_stats, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
         new=c("FS", "SOR", "MQRankSum", "ReadPosRankSum", "QD", "MQ", "DP"))

##remove rows with NA
variant_stats_nona <- na.omit(variant_stats)

variant_stats_nona_filtered <- subset(variant_stats_nona, filtration_status=="filtered")

########
## FS ## < 60??
########

ggplot(variant_stats_nona, aes(x=FS, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)+
  xlim(c(0, 100))

sum(variant_stats_nona_filtered$FS>60)

########
## QD ## > 2, NOT all above
########

ggplot(variant_stats_nona, aes(x=QD, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(variant_stats_nona_filtered$QD<2)

########
## DP ## > 10, none below
########
## Depth - remove lower than 10 and above 2*mean

ggplot(variant_stats_nona, aes(x=DP, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(variant_stats_nona_filtered$DP<10)

##mean of ???
mean(!is.na(unfSnps_variant_stats$DP))
##increases to 49.6 when less than 10 removed
dp_10 <- subset(variant_stats_nona, variant_stats_nona$DP>10)
mean(dp_10$DP)
##would reduce to 142,256
a <- subset(dp_10, dp_10$DP>100)
length(a$DP)

#########
## SOR ## GATK recc. > 3 ***
#########

ggplot(variant_stats_nona, aes(x=SOR, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(variant_stats_nona_filtered$SOR>3)

#########
## MQ ## GATK recc. > 40, all are
#########

ggplot(variant_stats_nona, aes(x=MQ, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(variant_stats_nona_filtered$MQ>40)/length(variant_stats_nona_filtered$MQ)

###############
## MQRankSum ## GATK recc. between -&+ 12.5 - all pass
###############

plot(density(variant_stats_nona$MQRankSum))

####################
## ReadPosRankSum ## GATK recc. between -&+ 8.0 - all pass
####################

ggplot(variant_stats_nona, aes(x=ReadPosRankSum, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

#### could still re-filter for depth>10 and SOR, all other filters pass