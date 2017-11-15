#!/bin/sh 

module load bcftools vcflib

# nb - you should run this from your own vcf directory

# randomly sampling the vcf
bcftools view $DATA_DIR/house_spanish_italian_chr15.vcf.gz | vcfrandomsample -r 0.1 > sparrow_subset.vcf

# only get biallelic data and fix counts
bcftools view -m2 -M2 sparrow_subset.vcf | vcffixup - > sparrow_subset_fix.vcf

# query the fixed data
bcftools query -f '%AC %AN %AF %QUAL %QD %DP %FS %MQ %BaseQRankSum %ReadPosRankSum\n' sparrow_subset_fix.vcf > snp_data.txt