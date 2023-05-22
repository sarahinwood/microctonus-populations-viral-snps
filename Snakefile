#!/usr/bin/env python3
import peppy

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

###########
# GLOBALS #
###########

#containers
gatk_container = 'docker://broadinstitute/gatk:4.3.0.0'
bcftools_container = 'docker://staphb/bcftools:1.16'
plink_container = 'docker://biocontainers/plink1.9:v1.90b6.6-181012-1-deb_cv1'

#########
# RULES #
#########

rule target:
    input:
       expand('output/01_gatk/{sample}_gatk_viral.g.vcf.gz', sample=all_samples),
       'output/01_gatk/MhV1_genotyped.vcf.gz',
       expand('output/02_filtering/unfiltered_{variant_type}s.vcf.gz', variant_type=["SNP", "INDEL"]),
       expand('output/02_filtering/{filtering}_SNPs_stats.txt', filtering=["unfiltered", "pass-filtered"]),
       'output/02_filtering/final_filtered_SNPs.vcf.gz',
       ##plink PCA
       expand('output/02_filtering/{location}_final_filtered_SNPs.vcf.gz', location=["SouthAmerican", "NZ"]),
       #expand('output/03_plink/no_ldpruning_{vcf_file}/mind_{mind}/filtered_snps_plink_pca.eigenvec', vcf_file=["final", "SouthAmerican_final", "NZ_final"], mind=["0.2", "0.5"]),
       #expand('output/03_plink/ld_pruned_{vcf_file}/mind_{mind}/filtered_snps_plink_pca.eigenvec', vcf_file=["final", "SouthAmerican_final", "NZ_final"], mind=["0.2", "0.5"])

###############
## plink PCA ##
###############

# pca on pruned
rule plink_pca_LD: # prunes down to 289 SNPs, and 133 samples by mind (69%) for all samples - what about historic?
    input:
        vcf = 'output/02_filtering/{vcf_file}_filtered_SNPs.vcf.gz',
        pruned = 'output/03_plink/ld_pruned_{vcf_file}/filtered_snps_plink.prune.in'
    output:
        'output/03_plink/ld_pruned_{vcf_file}/mind_{mind}/filtered_snps_plink_pca.eigenvec'
    params:
        out_dir = 'output/03_plink/ld_pruned_{vcf_file}/mind_{mind}/filtered_snps_plink_pca',
        mind = '{mind}'
    log:
        'output/logs/plink_pca_LD_{vcf_file}_{mind}.log'
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
        '--pca '
        '--mind  {paramns.mind} ' ## filter out samples with x% missing genotype data
        '--out {params.out_dir} '
        '&> {log}'

# prune dataset of variants in linkage - PCA relies on independent variables
## error as some samples present without genotype data
rule plink_prune_linkage:
    input:
        'output/02_filtering/{vcf_file}_filtered_SNPs.vcf.gz'
    output:
        'output/03_plink/ld_pruned_{vcf_file}/filtered_snps_plink.prune.in'
    params:
        indep = '10 10 0.1',     # 10 kb window, 10 bp step size, phased r2 threshold < 0.1 #### can definitely play around with this a bit more
        out = 'output/03_plink/ld_pruned_{vcf_file}/filtered_snps_plink'
    log:
        'output/logs/plink_prune_linkage_{vcf_file}.log'
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

##PCA on all snps - some likely in LD though
rule plink_pca_no_ldpruning: # down to 135 samples by mind (kept 70% samples) - which ones removed?
    input:
        vcf = 'output/02_filtering/{vcf_file}_filtered_SNPs.vcf.gz',
    output:
        'output/03_plink/no_ldpruning_{vcf_file}/mind_{mind}/filtered_snps_plink_pca.eigenvec'
    params:
        'output/03_plink/no_ldpruning_{vcf_file}/mind_{mind}/filtered_snps_plink_pca'
    log:
        'output/logs/plink_pca_no_ldpruning_{vcf_file}_{mind}.log'
    singularity:
        plink_container
    shell:
        '/usr/bin/plink1.9 '
        '--vcf {input.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--make-bed '
        '--pca '
        '--mind {params.mind} ' ## filter out samples with x% missing genotype data
        '--out {params} '
        '&> {log}'

#######################
## further filtering ##
#######################

## only historic SA samples - not historic hawkes bay
rule bcftools_split_timing:
    input:
        vcf = 'output/02_filtering/final_filtered_SNPs.vcf.gz',
        timing_samples = 'data/{location}_samples.txt'
    output:
        timing_vcf = 'output/02_filtering/{location}_final_filtered_SNPs.vcf.gz'
    log:
        "output/logs/bcftools_split_{location}.log"
    singularity:
        bcftools_container
    shell:
        'bcftools view -S {input.timing_samples} {input.vcf} -O z  -o {output.timing_vcf} 2> {log}'

# remove genotypes missing on more than X% of individuals and with minor allele freq. below Y
rule remove_singletons: #### no. SNPs retained: 8 ###
    input:
        'output/02_filtering/filtered_SNPs_ACAN.vcf.gz'
    output:
        'output/02_filtering/final_filtered_SNPs.vcf.gz' # quality filtered merged SNP file
    log:
        'output/logs/remove_singletons.log'
    singularity:
        bcftools_container
    shell:
        'bcftools filter '
        '-e "F_MISSING > 0.4 || MAF <= 0.05" '
        '-O z -o {output} {input} 2> {log}'
# when missing > 0.2 no. SNPs = 8
# when missing > 0.4 no. SNPs = 1238å

# AC==0 no alternative allele called
# AC==AN removes all sites where only alternative allele called - not true variant (assembly error)
# SnpGap removes variants close to indels - harder to call with certainty
# -m and -M set minimum and maximum number of alleles

# -m2 -M2 -v snps keeps only biallelic SNPs - minimum and maximum no. alleles is 2 aso only bi-allelic, indels removed ***** shouldn't be doing this for a haploid virus??
# q = minimum allele freq
rule filter_for_variants: #### no. SNPs retained: 2152 ###
    input:
        'output/02_filtering/pass-filtered_SNPs.vcf.gz'
    output:
        'output/02_filtering/filtered_SNPs_ACAN.vcf.gz'
    log:
        'output/logs/filter_for_variants.log'
    singularity:
        bcftools_container
    shell:
        'bcftools filter -e "AC==0 || AC==AN" --SnpGap 10 {input} | bcftools view -q 0.1 -m2 -M10 -v snps -O z -o {output} 2> {log}' # max should also be 2, not pooled anymore

######################
## filtration stats ##
######################

rule variant_stats:
    input:
        'output/02_filtering/{filtering}_SNPs.vcf.gz'
    output:
        'output/02_filtering/{filtering}_SNPs_stats.txt'
    log:
        'output/logs/{filtering}_variant_stats.log'
    singularity:
        bcftools_container
    shell:
        'bcftools query {input} -f "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\t%AC\n" > {output} 2> {log}'

#############################
## gatk variant filtration ##
#############################

# then select only variants that passed - i.e. didn't meet filters below
rule select_pass_merged: #### no. SNPs retained: 3541 ###
    input:
        'output/02_filtering/filtered_SNPs.vcf.gz'
    output:
        'output/02_filtering/pass-filtered_SNPs.vcf.gz'
    singularity:
        bcftools_container
    shell:
        'bcftools view -f "PASS" {input} -o {output} -Oz'

# filter with GATK recc.s - https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# filters will get rid of any SNPs that meet the below rules
rule filter_SNPs:
    input:
        'output/02_filtering/unfiltered_SNPs.vcf.gz'
    output:
        'output/02_filtering/filtered_SNPs.vcf.gz'
    log:
        'output/logs/filter_SNPs.log'
    singularity:
        gatk_container
    shell:
        'gatk VariantFiltration '
        '-V {input} '
        '-filter "QUAL < 30.0" --filter-name "QUAL30" ' # recc.d in some GATK articles but not all
        '-filter-expression "DP < 50.0" --filter-name "DP50.0"  ' # filter by depth - filters by depth across all samples not by each sample individually, can probably be stricter here
        '-filter "FS > 60.0" --filter-name "FS60" ' # fisher score for strand bias
        '-filter "SOR > 3.0" --filter-name "SOR3" ' # strand odds ratio for strand bias
        '-filter "MQ < 40.0" --filter-name "MQ40" ' # root mean square mapping quality
        '-filter "MQRankSum < -12.5" --filter-name "MQRankSum12.5" ' #rank sum test for mapping quality
        '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum8" ' # rank sum test for site position
        '-filter "AC < 2 " --filter-name "AC2" ' # minimum allele count filter - at least 2 - shouldn't be using this - should always = 2 for diploid, but for haploid virus which might have multiple variants in an individual shouldn't use
        '--missing-values-evaluate-as-failing true '
        '-O {output} '
        '2> {log}'

# Tristan used:
    # min-mac 2 - minimum minor allele count, also min/max allele count = 2, I've done this later   ***** need to change max to 2 as well - have it incorrectly later on from pooled genome sample - or should these actually be 1 because not diploid
    # minDP of 200, max of 10k                                                                      ***** basically no difference, use Tristan's filters
    # QD > 10, I've not used this!!!!!                                                              ***** can trial both 2 and 10
    # MQ > 50, I've used 40                                                                         ***** this increase makes a big difference, trial 40 and 50
    # no indels or multiallelics                                                                    ***** then my extra remaining filters don't change too much so probably keep them
    # minimum allele freq = 0.1:minor - I've done this later                                        ***** need to think about minimum allele count filter to ensure it makes sense

rule select_variants: #### no. SNPs: 6855 ###
    input:
        'output/01_gatk/MhV1_genotyped.vcf.gz'
    output:
        'output/02_filtering/unfiltered_{variant_type}s.vcf.gz'
    params:
        vartype = '{variant_type}'
    log:
        'output/logs/select_only_{variant_type}s.log'
    singularity:
        gatk_container
    shell:
        'gatk SelectVariants '
        '-V {input} '
        '-select-type {params.vartype} '
        '-O {output} '
        '2> {log}'

##########################
## MhFV variant calling ##
##########################

rule genotype_gvcfs_combined: #### no. variants: 9563 ###
    input:
        genome = 'data/Mh_MhV1.fa',
        vcf = 'output/01_gatk/MhV1_combined.g.vcf.gz'
    output:
        'output/01_gatk/MhV1_genotyped.vcf.gz'
    log:
        'output/logs/genotype_gvcfs.log'
    singularity:
        gatk_container
    threads:
        40
    shell:
        ' gatk --java-options "-Xmx4g" GenotypeGVCFs '
        '-R {input.genome} '
        '-V {input.vcf} '
        '-O {output} '
        '2> {log}'

rule combine_gvcfs:
    input:
        genome = 'data/Mh_MhV1.fa',
        viral = expand('output/01_gatk/{sample}_gatk_viral.g.vcf.gz', sample=all_samples)
    output:
        combined = 'output/01_gatk/MhV1_combined.g.vcf.gz'
    params:
        wd = 'output/01_gatk/',
        vcf_list = 'output/01_gatk/gVCF.list'
    log:
        'output/logs/combine_gvcfs.log'
    singularity:
        gatk_container
    threads:
        40
    shell:
        'find {params.wd} -type f -name "*_gatk_viral.g.vcf.gz" > {params.vcf_list} || exit 1 ; '
        'gatk CombineGVCFs '
        '-R {input.genome} '
        '--variant {params.vcf_list} '
        '-O {output.combined} '
        '2> {log}'

rule gatk_viral:
    input:
        genome = 'data/Mh_MhV1.fa',
        genome_dict = 'data/Mh_MhV1.dict',
        genome_index = 'data/Mh_MhV1.fa.fai',
        bam = 'output/01_gatk/{sample}_rgs_dd_s.bam',
        index = 'output/01_gatk/{sample}_rgs_dd_s.bam.bai'
    output:
        gvcf = 'output/01_gatk/{sample}_gatk_viral.g.vcf.gz'
    threads:
        40
    log:
        'output/logs/gatk_viral_{sample}.log'
    singularity:
        gatk_container
    shell:
        'gatk --java-options "-Xmx4g" HaplotypeCaller '
        '-R {input.genome} '
        '-I {input.bam} '
        '-L scaffold_1_MhV1 ' # restrict to only viral contig
        '--ploidy 1 ' # haploid for virus
        '--native-pair-hmm-threads {threads} '
        '-ERC GVCF ' # output gVCF files
        '-O {output.gvcf} '
        '2> {log}'

# Tristan never did Base Quality Score Recalibration (BQSR) but it is reccommended on GATK website so may be worth a try

rule samtools_index_gatk:
    input:
        bam = 'output/01_gatk/{sample}_rgs_dd_s.bam'
    output:
        index = 'output/01_gatk/{sample}_rgs_dd_s.bam.bai'
    log:
        'output/logs/samtools_index_gatk_{sample}.log'
    threads:
        5
    shell:
        'samtools index '
        '{input.bam} '
        '2> {log}'

rule gatk_sortsam:
    input:
        'output/01_gatk/{sample}_rgs_dd.bam'
    output:
        'output/01_gatk/{sample}_rgs_dd_s.bam'
    log:
        'output/logs/gatk_sortsam_{sample}.log'
    singularity:
        gatk_container
    threads:
        5
    shell:
        'gatk --java-options "-Xmx8g" SortSam '
        'I={input} '
        'O={output} '
        'SORT_ORDER=coordinate '
        '2> {log}'

rule gatk_mark_dups:
    input:
        bam = 'output/01_gatk/{sample}_rgs.bam'
    output:
        md_bam = temp('output/01_gatk/{sample}_rgs_dd.bam'),
        metrics = 'output/01_gatk/{sample}_mark_dups_metrics.txt'
    log:
        'output/logs/gatk_mark_dups_{sample}.log'
    singularity:
        gatk_container
    threads:
        5
    shell:
        'gatk --java-options "-Xmx8g" MarkDuplicates '
        'I={input.bam} '
        'ASSUME_SORT_ORDER=coordinate '
        'O={output.md_bam} '
        'M={output.metrics} '
        '2> {log}'

rule gatk_readgroups:
    input:
        'output/samtools/{sample}_sorted_fixed.bam'
    output:
        temp('output/01_gatk/{sample}_rgs.bam')
    params:
        sample_id = '{sample}'
    log:
        'output/logs/gatk_readgroups_{sample}.log'
    singularity:
        gatk_container
    threads:
        5
    shell:
        'gatk --java-options "-Xmx8g" AddOrReplaceReadGroups '
        'I={input} O={output} '
        'RGID={params.sample_id} RGLB={params.sample_id} RGPU={params.sample_id} RGSM={params.sample_id} RGPL=ILLUMINA '
        '2> {log}'

#######################
## fix bams for GATK ##
#######################

### reads with 0x2 flag cause errors in MarkDuplicates ###
rule samtools_fix_bams:
    input:
        'output/samtools/{sample}_sorted.bam'
    output:
        temp('output/samtools/{sample}_sorted_fixed.bam')
    threads:
        5
    shell:
        'samtools view -f 0x2 {input} -b > {output}'

rule sort_bams_coordinate:
    input:
        'output/samtools/{sample}_fixed.bam'
    output:
        temp('output/samtools/{sample}_sorted.bam')
    threads:
        5
    shell:
        'samtools sort {input} -o {output}'

rule fix_read_mates:
    input:
        'output/samtools/{sample}_name_sorted.bam'
    output:
        temp('output/samtools/{sample}_fixed.bam')
    threads:
        5
    shell:
        'samtools fixmate {input} {output}'

rule sort_bams_name:
    input:
        'output/bwa_Mh_MhFV/{sample}_sorted.bam'
    output:
        temp('output/samtools/{sample}_name_sorted.bam')
    threads:
        5
    shell:
        'samtools sort -n {input} -o {output}'

#########################
## prep fasta for GATK ##
#########################

rule gatk_seq_dict:
    input:
        'data/Mh_MhV1.fa'
    output:
        'data/Mh_MhV1.dict'
    log:
        'output/logs/gatk_seq_dict.log'
    singularity:
        gatk_container
    shell:
        'gatk CreateSequenceDictionary -R {input} 2> {log}'

rule faidx:
    input:
        fasta = 'data/Mh_MhV1.fa'
    output:
        fai = 'data/Mh_MhV1.fa.fai'
    log:
        'output/logs/faidx.log'
    shell:
        'samtools faidx {input.fasta} -o {output.fai} 2> {log} '
