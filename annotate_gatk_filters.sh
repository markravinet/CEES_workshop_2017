#!/bin/bash

module load bcftools vcftools

## define variables
PICARD_PATH=/projects/cees/bin/picard/2.7.1/picard.jar
GATK_PATH=/cluster/software/VERSIONS/gatk-3.7/GenomeAnalysisTK.jar
REF=/work/users/msravine/ref/house_sparrow_genome_assembly-18-11-14_masked.fa

INPUT_VCF=house_spanish_italian_chr15.vcf.gz
SNP_VCF=house_spanish_italian_chr15_snps_only.vcf
ANNOTATE_VCF=house_spanish_italian_chr15_snps_only_ann.vcf

# first things first, we want to remove indels (to make life easier)
# also this creates an uncompressed vcf - which GATK requires
vcftools --gzvcf $INPUT_VCF --remove-indels --recode --recode-INFO-all \
 --stdout > $SNP_VCF
 
# now we use GATK to annotate the vcf
echo "Annotating SNPs"
java -Xmx2g -jar $GATK_PATH -T VariantFiltration \
-R $REF \
-V $SNP_VCF \
--filterExpression "QD < 5.0 || FS > 10.0 || MQ < 40.0 || MQRankSum < -2.5 || ReadPosRankSum < -2.0" \
--filterName "hard_pass_filter" \
-o $ANNOTATE_VCF

# bgzip and index
bgzip -f $ANNOTATE_VCF

bcftools index ${ANNOTATE_VCF}.gz


