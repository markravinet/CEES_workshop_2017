#!/bin/sh

module load vcftools bcftools

cd /work/users/msravine/course_prep

VCF=/projects/cees2/in_progress/cees_workshop/vcf/house_spanish_italian_chr15.vcf.gz
# VCF=/projects/cees2/in_progress/for_workshop/angsd/HC_GP_GL_vcf_v4.2_0.95.vcf.gz
#VCF=/work/users/olekto/sparrow_course/house_sparrow_genome_assembly-18-11-14_masked_seg_15/house_sparrow_genome_assembly-18-11-14_masked_seg_15_sparrow_vars_freebayes.vcf.gz

# mean depth per ind
# mean depth per variant
# heterozygosity
# major allele freq
# site missingness
# ind missingness

# get site identities
# bcftools view -H $VCF | cut -f1-2 > chr15_snps

# allele freq
vcftools --gzvcf $VCF --remove-indels --min-alleles 2 --max-alleles 2 --freq2 --out chr15

# depth per ind
vcftools --gzvcf $VCF --remove-indels --min-alleles 2 --max-alleles 2 --depth --out chr15

# depth per site
vcftools --gzvcf $VCF --remove-indels --min-alleles 2 --max-alleles 2F --site-mean-depth --out chr15

# het per ind
vcftools --gzvcf $VCF --remove-indels --min-alleles 2 --max-alleles 2 --het --out chr15

# site quality
vcftools --gzvcf $VCF --remove-indels --min-alleles 2 --max-alleles 2 --site-quality --out chr15

# missing ind
vcftools --gzvcf $VCF --remove-indels --min-alleles 2 --max-alleles 2 --missing-indv --out chr15

# missing site
vcftools --gzvcf $VCF --remove-indels --min-alleles 2 --max-alleles 2 --missing-site --out chr15