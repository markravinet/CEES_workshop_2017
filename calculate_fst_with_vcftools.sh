#!/bin/sh

module load vcftools bcftools

VCF=/projects/cees2/in_progress/cees_workshop/vcf/house_spanish_italian_chr15_filtered.vcf.gz
WINDOW=100000
STEP=25000

# generate pop files
bcftools query -l $VCF | grep "house" > house
bcftools query -l $VCF | grep "spanish" > spanish

vcftools --gzvcf $VCF \
--fst-window-size $WINDOW --fst-window-step $STEP \
--weir-fst-pop house --weir-fst-pop spanish --out house_spanish
