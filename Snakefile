#!/usr/bin/env python3
import peppy
import pathlib2

###########
# GLOBALS #
###########

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']


## sample subsets ##

# Function to load samples from a text file
def load_samples(filename):
    with open(filename) as f:
        return [line.strip() for line in f if line.strip()]

# Load the sample list from the text file
all_hyp = load_samples("data/sample_lists/all_hyp.txt")
Mhyp_75cov = load_samples("data/sample_lists/Mhyp_75cov.txt")

## vcf subsets and lists
vcf_subsets = ["Mhyp_75cov_only"] #, "Mhyp_85cov_only", "all", "Mhyp_only"
location_splits =["SouthAmerican", "NZ", "NZ_withtom", "NZ_postrelease_contemp", "NZ_postrelease_early", "Mhyp", "SouthAmericanCohort4", "SouthAmericanCohort4_earlyCohort5"] #"Chatham_PCA_75cov"

###########
# GLOBALS #
###########

#containers
gatk_container = 'docker://broadinstitute/gatk:4.3.0.0'
picard_container = 'docker://broadinstitute/picard:2.27.5'
bcftools_container = 'docker://staphb/bcftools:1.16'
plink_container = 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'
plink_v2_container = 'library://sinwood/plink/plink_2.0:0.0.1'
admixture_container = 'docker://evolbioinfo/admixture:v1.3.0'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.19:0.0.2'
samtools_container = 'docker://staphb/samtools:1.10'

#########
# RULES #
#########

rule target:
    input:
            ## variant calling and genotyping
       #expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_unfiltered_{subset_vcf}_{variant_type}s.vcf.gz', variant_type=["SNP", "INDEL"], subset_vcf=vcf_subsets),
       #expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/unfiltered_SNPs_stats_{subset_vcf}.txt', subset_vcf=vcf_subsets),
            ## filtering
       #expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/filtered_SNPs_{subset_vcf}.vcf.gz', subset_vcf=vcf_subsets),
       #expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/gatk_pass-filtered_SNPs_{subset_vcf}.vcf.gz', subset_vcf=vcf_subsets),
            ## filtering stats
       #expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/stats/{subset_vcf}_SNPs_stats.txt', subset_vcf=vcf_subsets),
       #expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/stats/bcftools_stats_filtered_SNPs_{subset_vcf}.txt', subset_vcf=vcf_subsets),
            ## split VCFs
       expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/04_filtered_split_vcfs/{location}_Mhyp_75cov_only_filtered_SNPs.vcf.gz', location=location_splits),
            ##plink PCA
        expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/no_ldpruning_{vcf_file}_Mhyp_75cov_only/mind_{mind}/filtered_SNPs_plink_pca.eigenvec', vcf_file=["SouthAmerican", "SouthAmericanCohort4", "NZ", "Mhyp", "SouthAmericanCohort4_earlyCohort5"], mind=["0.2"]),
       expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{vcf_file}_Mhyp_75cov_only/mind_{mind}/filtered_SNPs_plink_pca.eigenvec', vcf_file=["SouthAmerican", "SouthAmericanCohort4", "SouthAmericanCohort4_earlyCohort5", "NZ", "Mhyp"], mind=["0.2"]),
       #      ## population-specific SNP ID
       #  expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/population_specific_SNPs/ecotype/ecotype.{ecotype}.afreq', ecotype=["Ascasubi", "RioNegro", "Mendoza"]),
       #  expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/population_specific_SNPs/line/line.{line}.afreq', line=["AS28", "AS6", "MN1", "RN36"]),
       #  expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/population_specific_SNPs/NZ_postrelease_early/location.{location}.afreq', location=["HawkesBay", "Reporoa", "Ruakura", "Wellsford"]),
       #  expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/population_specific_SNPs/NZ_postrelease_contemp/location.{location}.afreq', location=["ChathamIsland", "Featherston", "Fortrose", "Lincoln", "Mararoa", "Mossburn", "Ruakaka", "Ruakura", "Timaru", "WestCoastSite2"]),
       #      # proportions of ecotype-specific SNPs
       #  expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/population_specific_SNPs/ecotype/sample_proportions/{sample}/{c4_location}/README.txt', sample='Ascasubi_Argentina_1997_1', c4_location='Ascasubi')

###############
## plink PCA ##
###############

# pca on pruned
rule plink_pca_LD: # prunes down to ?? SNPs, and ?? samples by mind (69%) for all samples - what about historic?
    input:
        vcf = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/04_filtered_split_vcfs/{location}_{subset_vcf}_filtered_SNPs.vcf.gz',
        pruned = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/filtered_SNPs_plink.prune.in'
    output:
        eigenvec = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/mind_{mind}/filtered_SNPs_plink_pca.eigenvec',
        bed = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/mind_{mind}/filtered_SNPs_plink_pca.bed',
        bim = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/mind_{mind}/filtered_SNPs_plink_pca.bim',
        fam = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/mind_{mind}/filtered_SNPs_plink_pca.fam'
    params:
        out_dir = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/mind_{mind}/filtered_SNPs_plink_pca',
        mind = '{mind}'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/plink_pca_LD_{location}_{mind}_{subset_vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--vcf {input.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--extract {input.pruned} '
        '--make-bed '
        '--pca var-wts '
        '--mind  {params.mind} ' ## filter out samples with x% missing genotype data
        '--out {params.out_dir} '
        '&> {log}'

# prune dataset of variants in linkage - PCA relies on independent variables
## error as some samples present without genotype data
rule plink_prune_linkage:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/04_filtered_split_vcfs/{location}_{subset_vcf}_filtered_SNPs.vcf.gz'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/filtered_SNPs_plink.prune.in'
    params:
        indep = '10 10 0.1',     # 10 kb window, 10 bp step size, phased r2 threshold < 0.1 #### can definitely play around with this a bit more
        out = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/ld_pruned_{location}_{subset_vcf}/filtered_SNPs_plink'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/plink_prune_linkage_{location}_{subset_vcf}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--vcf {input} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--indep-pairwise {params.indep} '
        '--out {params.out} '
        '&> {log}'

##PCA on all SNPs - some likely in LD though
rule plink_pca_no_ldpruning: # down to 135 samples by mind (kept 70% samples) - which ones removed?
    input:
        vcf = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/04_filtered_split_vcfs/{location}_{subset_vcf}_filtered_SNPs.vcf.gz',
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/no_ldpruning_{location}_{subset_vcf}/mind_{mind}/filtered_SNPs_plink_pca.eigenvec'
    params:
        out_dir = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/05_plink/no_ldpruning_{location}_{subset_vcf}/mind_{mind}/filtered_SNPs_plink_pca',
        mind = '{mind}'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/plink_pca_no_ldpruning_{location}_{subset_vcf}_{mind}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--vcf {input.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--make-bed '
        '--pca var-wts '
        '--mind {params.mind} ' ## filter out samples with x% missing genotype data
        '--out {params.out_dir} '
        '&> {log}'

## what effect is mind having here when only feeding in samples with >75% coverage?

#########################
## filter & split VCFs ##
#########################

rule bcftools_split_timing:
    input:
        vcf = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/filtered_SNPs_Mhyp_75cov_only.vcf.gz',
        timing_samples = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/sample_lists/{location}_samples.txt'
    output:
        timing_vcf = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/04_filtered_split_vcfs/{location}_Mhyp_75cov_only_filtered_SNPs.vcf.gz'
    log:
        "/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/bcftools_split_{location}_Mhyp_75cov_only.log"
    singularity:
        bcftools_container
    shell:
        'nice bcftools view -S {input.timing_samples} {input.vcf} --force-samples -O z  -o {output.timing_vcf} 2> {log}'

# stats on filtered snps
rule bcf_stats_filtered:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/filtered_SNPs_{subset_vcf}.vcf.gz'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/stats/bcftools_stats_filtered_SNPs_{subset_vcf}.txt'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/bcf_stats_filtered_{subset_vcf}.log'
    singularity:
        bcftools_container  
    shell:
        'nice bcftools stats {input} '
        '--samples '
        '"-" '
        '> {output} '
        '2> {log}'

rule variant_stats_filtered:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/filtered_SNPs_{subset_vcf}.vcf.gz'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/stats/{subset_vcf}_SNPs_stats.txt'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/{subset_vcf}_variant_stats.log'
    singularity:
        bcftools_container
    shell:
        'bcftools query {input} -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\t%AC\n" > {output} 2> {log}'

# filtering based on Josephs suggestions #### no. SNPs retained: ?? was previously 3541 ###
rule bcftools_filtering:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_unfiltered_{subset_vcf}_SNPs.vcf.gz'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/03_filtering/filtered_SNPs_{subset_vcf}.vcf.gz'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/bcf_filtering_vcf_{subset_vcf}.log'
    singularity:
        bcftools_container
    shell:
        'nice bcftools view {input} '
        '-i "AVG(GQ)>20 & QUAL>40 & F_MISSING<0.2" '    # include SNPs where F_MISSING <0.2 (could be to 0.3, 0.4), quality>40, & average of genotype qual.s larger than 20 (GQ=20 is 1% chance call is incorrect)
        '-q 0.05 '                                      # min. allele freq, could do 0.02, 0.01
        '-O z -o {output} '                             # output file name
        '2> {log}'
## other things I could but don't know if I should filter on:
    # could try qual>30 and see what difference that makes - could be too stringent for a virus
    # AC==0 - remove any variants with allele count of zero (alt allele not present)
    # AC==AN - remove any variants where allele count is equal to allele number - not a true SNP, is likely an assembly error
    # depth filters - min and max - what influence does this make when feeding in only samples with >75% coverage? Tristan had used min 200 max 10k, have a read/think/play around

#################################################
## separate SNPs & INDELs --> unfiltered stats ##
#################################################

rule variant_stats_unfiltered:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_unfiltered_{subset_vcf}_SNPs.vcf.gz'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/unfiltered_SNPs_stats_{subset_vcf}.txt'
    singularity:
        bcftools_container
    shell:
        'bcftools query {input} -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\t%AC\t[%GQ]\n" > {output}'

## separate SNPs and INDELs
rule select_variants: #### no. SNPs: 7537 (was 6855) ###
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_unfiltered_genotyped_{subset_vcf}.vcf.gz'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_unfiltered_{subset_vcf}_{variant_type}s.vcf.gz'
    params:
        vartype = '{variant_type}'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/select_only_{subset_vcf}_{variant_type}s.log'
    singularity:
        gatk_container
    shell:
        'nice gatk SelectVariants '
        '-V {input} '
        '-select-type {params.vartype} '
        '-O {output} '
        '2> {log}'

##############################################
## Genotype GVCFs -  all & >75% cov samples ##
##############################################

## all samples and then only high coverage samples to genotype separately - low coverage samples likely to add noise to analysis
rule genotype_gvcfs_combined: #### no. variants: 10699 was 9563 ###
    input:
        genome = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa',
        vcf = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_combined_{subset_vcf}.g.vcf.gz'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_unfiltered_genotyped_{subset_vcf}.vcf.gz'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/genotype_gvcfs_{subset_vcf}.log'
    singularity:
        gatk_container
    shell:
        'nice gatk --java-options "-Xmx12g" GenotypeGVCFs '
        '-R {input.genome} '
        '-V {input.vcf} '
        '-O {output} '
        '2> {log}'

rule combine_gvcfs_Mhyp_75cov_only:
    input:
        genome = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa',
        viral = expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/01_gatk/{sample}_gatk_viral.g.vcf.gz', sample=Mhyp_75cov),
        vcf_list = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/sample_lists/Mhyp_75cov_vcfs.txt'
    output:
        combined = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_combined_Mhyp_75cov_only.g.vcf.gz'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/combine_gvcfs_Mhyp_75cov.log'
    singularity:
        gatk_container
    shell:
        'nice gatk CombineGVCFs '
        '-R {input.genome} '
        '--variant {input.vcf_list} '
        '-O {output.combined} '
        '2> {log}'

rule combine_gvcfs_all_hyp:
    input:
        genome = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa',
        viral = expand('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/01_gatk/{sample}_gatk_viral.g.vcf.gz', sample=all_hyp),
        vcf_list = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/sample_lists/all_hyp_vcfs.txt'
    output:
        combined = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/02_gatk_combined/MhFV_combined_Mhyp_only.g.vcf.gz'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/combine_gvcfs_all_hyp.log'
    singularity:
        gatk_container
    shell:
        'nice gatk CombineGVCFs '
        '-R {input.genome} '
        '--variant {input.vcf_list} '
        '-O {output.combined} '
        '2> {log}'

######################################
## haplotype caller for each sample ##
######################################

rule gatk__haplotypecaller_viral:
    input:
        genome = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa',
        genome_dict = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.dict',
        genome_index = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa.fai',
        bam = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_markdup.cram',
        index = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_markdup.cram.crai'
    output:
        gvcf = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/01_gatk/{sample}_gatk_viral.g.vcf.gz'
    threads:
        5
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/gatk_viral_{sample}.log'
    singularity:
        gatk_container
    shell:
        'nice gatk --java-options "-Xmx12g" HaplotypeCaller '
        '-R {input.genome} '
        '-I {input.bam} '
        '-L MhFV_scaffold ' # restrict to only viral contig
        '--ploidy 1 ' # haploid for virus
        '--native-pair-hmm-threads {threads} '
        '-ERC GVCF ' # output gVCF files
        '-O {output.gvcf} '
        '2> {log}'

#################################
## transform bams and markdups ##
#################################

rule samtools_index_gatk:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_markdup.cram'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_markdup.cram.crai'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/samtools_index_gatk_{sample}.log'
    singularity:
        samtools_container
    shell:
        'nice samtools index '
        '{input} '
        '2> {log}'

rule samtools_markdup:
    input:
        cram = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_fixed_sorted.cram',
        fa = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_markdup.cram'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/samtools_markdup/{sample}.log'
    threads:
        2
    singularity:
        samtools_container
    shell:
        'nice samtools markdup '
        '--threads {threads} '
        '--reference {input.fa} '
        '-O CRAM {input.cram} '
        '{output} '
        '2> {log}'

# markdups needs position order
rule sort_bams_position:
    input:
        cram = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_fixed.cram',
        fa = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa'
    output:
        temp('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_fixed_sorted.cram')
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/sort_bams_position/{sample}.log'
    threads:
        2
    singularity:
        samtools_container
    shell:
        'nice samtools sort '
        '--threads {threads} '
        '--reference {input.fa} '
        '{input.cram} '
        '-o {output} '
        '2> {log}'

# based off Josephs nf
# fills in mate coordinates, -r removing secondary and unmapped reads, -c adding template cigar ct tag, -m adds mate score tag used by markdup
rule fix_read_mates:
    input:
        cram = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_name_sorted_cr.cram',
        fa = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa'
    output:
        temp('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_fixed.cram')
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/fix_read_mates/{sample}.log'
    threads:
        2
    singularity:
        samtools_container
    shell:
        'nice samtools fixmate '
        '-r -c -m '
        '--threads {threads} '
        '--reference {input.fa} '
        '{input.cram} '
        '{output} '
        '2> {log}'

# sorting by name
rule sort_bams_name:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_sorted_cr.cram'
    output:
        temp('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_name_sorted_cr.cram')
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/sort_bams_name/{sample}.log'
    threads:
        2
    singularity:
        samtools_container
    shell:
        'nice samtools sort '
        '-n -O CRAM '
        '--threads {threads} '
        '{input} '
        '-o {output} '
        '2> {log}'

# converting to CRAM to help speed up following steps
rule convert_bamcram:
    input:
        bam = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_sorted_rgs.bam',
        fa = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa',
        fai = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa.fai'
    output:
        temp('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_sorted_cr.cram')
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/convert_bamcram/{sample}.log'
    threads:
        2
    container:
        samtools_container
    shell:
        'nice samtools view '
        '-C ' #cram output
        '--reference {input.fa} '
        '-t {input.fai} '
        '--threads {threads} '
        '{input.bam} '
        '-o {output} '
        '2> {log}'

rule gatk_readgroups:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam'
    output:
        temp('/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/00_crams/{sample}_sorted_rgs.bam')
    params:
        sample_id = '{sample}'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/gatk_readgroups/{sample}.log'
    singularity:
        gatk_container
    shell:
        'nice gatk --java-options "-Xmx8g" AddOrReplaceReadGroups '
        'I={input} O={output} '
        'RGID={params.sample_id} RGLB={params.sample_id} RGPU={params.sample_id} RGSM={params.sample_id} RGPL=ILLUMINA '
        '2> {log}'

#############################
## prep fasta ref for GATK ##
#############################

rule gatk_seq_dict:
    input:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa'
    output:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.dict'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/gatk_seq_dict.log'
    singularity:
        gatk_container
    shell:
        'gatk CreateSequenceDictionary -R {input} 2> {log}'

rule faidx:
    input:
        fasta = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa'
    output:
        fai = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/data/Mh_MhmtDNA_MhFV.fa.fai'
    log:
        '/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-populations-viral-snps/output/logs/faidx.log'
    shell:
        'samtools faidx {input.fasta} -o {output.fai} 2> {log} '

