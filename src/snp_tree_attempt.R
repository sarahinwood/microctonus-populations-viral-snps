library(SeqArray)
library(SNPRelate)

gvcf_path <- "output/gatk/MhV1_genotyped.vcf.gz"
gds_path <- "output/gatk/MhV1_genotyped.gds"

seqVCF2GDS(gvcf_path, gds_path,
           ignore.chr.prefix = "scaffold_1_MhV", parallel=2)
# info.import = "", fmt.import = "", genotype.var.name = "",

snpgdsSummary(snpgdsExampleFileName())
genofile <- snpgdsOpen(snpgdsExampleFileName())
get.attr.gdsn(index.gdsn(genofile, "snp.chromosome"))
snpgdsSummary("test.gds")