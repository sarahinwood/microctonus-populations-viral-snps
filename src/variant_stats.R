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
#filSnps_variant_stats <- fread("output/02_filtering/pass-filtered_SNPs_stats.txt", header=F, na.strings=".")
#filSnps_variant_stats$filtration_status <- paste("filtered")
#join tables
#variant_stats <- full_join(unfSnps_variant_stats, filSnps_variant_stats)

unfSnps_variant_stats <- fread("output/02_filtering/unfiltered_SNPs_stats.txt", header=F, na.strings=".")
unfSnps_variant_stats$filtration_status <- paste("unfiltered")
setnames(unfSnps_variant_stats, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("CHROM", "POS", "REF", "ALT", "QUAL", "FS", "SOR", "MQRankSum", "ReadPosRankSum", "QD", "MQ", "DP", "AC"))
##remove rows with NA
unfvariant_stats_nona <- na.omit(unfSnps_variant_stats)

# DP - barely any difference here - can change to Tristans
DP_Tristan_filters <- subset(unfvariant_stats_nona, DP>200 & DP<10000)
DP_My_filters <- subset(unfvariant_stats_nona, DP>50)
paste("Tristan's filter removes an extra", length(DP_My_filters$FS)-length(DP_Tristan_filters$FS), "SNPs", sep=" ")

# QD - moderate difference (442 SNPs) not currently implemented but should be, trial both?
QD_Tristan_filters <- subset(unfvariant_stats_nona, QD>10)
QD_My_filters <- subset(unfvariant_stats_nona, QD>2)
paste("Tristan's filter removes an extra", length(QD_My_filters$FS)-length(QD_Tristan_filters$FS), "SNPs", sep=" ")

# MQ - makes sizeable difference (1030 SNPs) and should trial both?
MQ_Tristan_filters <- subset(unfvariant_stats_nona, MQ>50)
MQ_My_filters <- subset(unfvariant_stats_nona, MQ>40)
paste("Tristan's filter removes an extra", length(MQ_My_filters$FS)-length(MQ_Tristan_filters$FS), "SNPs", sep=" ")

Tristan_filters <- subset(unfvariant_stats_nona, DP>200 & DP<10000 & QD>10 & MQ>50)

My_filters <- subset(unfvariant_stats_nona, QUAL>30 & DP>50 & FS<60 & SOR<3 & MQ>40 & MQRankSum>-12.5 & ReadPosRankSum>-8.0 & AC>=2)
My_filters_QD <- subset(unfvariant_stats_nona, QD>2 & DP>50 & FS<60 & SOR<3 & MQ>40 & MQRankSum>-12.5 & ReadPosRankSum>-8.0)
My_filters_QD_10 <- subset(unfvariant_stats_nona, QD>10 & DP>50 & FS<60 & SOR<3 & MQ>40 & MQRankSum>-12.5 & ReadPosRankSum>-8.0)

Tristan_and_my_other_filters <- subset(unfvariant_stats_nona, QD>10 & DP>200 & DP<10000 & FS<60 & SOR<3 & MQ>50 & MQRankSum>-12.5 & ReadPosRankSum>-8.0)
# adding my extra filters doesn't make a large difference so may as well

########
## QD ## > 2, NOT all above
########

ggplot(unfvariant_stats_nona, aes(x=QD, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(unfvariant_stats_nona_filtered$QD<2)

########
## FS ## < 60??
########

ggplot(unfvariant_stats_nona, aes(x=FS, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)+
  xlim(c(0, 100))

sum(unfvariant_stats_nona_filtered$FS>60)

########
## DP ## > 10, none below
########
## Depth - remove lower than 10 and above 2*mean

ggplot(unfvariant_stats_nona, aes(x=DP, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(unfvariant_stats_nona_filtered$DP<10)

##mean of ???
mean(!is.na(unfSnps_variant_stats$DP))
##increases to 49.6 when less than 10 removed
dp_10 <- subset(unfvariant_stats_nona, variant_stats_nona$DP>10)
mean(dp_10$DP)
##would reduce to 142,256
a <- subset(dp_10, dp_10$DP>100)
length(a$DP)

#########
## SOR ## GATK recc. > 3 ***
#########

ggplot(unfvariant_stats_nona, aes(x=SOR, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(unfvariant_stats_nona_filtered$SOR>3)

#########
## MQ ## GATK recc. > 40, all are
#########

ggplot(unfvariant_stats_nona, aes(x=MQ, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

sum(unfvariant_stats_nona_filtered$MQ>40)/length(unfvariant_stats_nona_filtered$MQ)

###############
## MQRankSum ## GATK recc. between -&+ 12.5 - all pass
###############

plot(density(unfvariant_stats_nona$MQRankSum))

####################
## ReadPosRankSum ## GATK recc. between -&+ 8.0 - all pass
####################

ggplot(unfvariant_stats_nona, aes(x=ReadPosRankSum, fill=filtration_status))+
  geom_density()+
  facet_wrap(~filtration_status)

#### could still re-filter for depth>10 and SOR, all other filters pass