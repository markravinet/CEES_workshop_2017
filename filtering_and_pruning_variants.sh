#!/bin/sh

module load bcftools plink/1.90b3b vcftools

VCF=/projects/cees2/in_progress/cees_workshop/vcf/house_spanish_italian_chr15_snps_only_ann.vcf.gz
FILTERED_VCF=/projects/cees2/in_progress/cees_workshop/vcf/house_spanish_italian_chr15_filtered.vcf.gz

# throwout sites with spanning deletions
bcftools view -e 'ALT="*"' -O z -o ${VCF%.*.*}_sd.vcf.gz $VCF

# redefine the variable
VCF=${VCF%.*.*}_sd.vcf.gz

# first we will have a quick look at the annotated vcf
bcftools view -H $VCF | head

# next we need to filter the vcf
vcftools --gzvcf $VCF --remove-filtered-all --remove-indels \
--min-alleles 2 --max-alleles 2 \
--max-missing 0.7 --min-meanDP 10 \
--recode --recode-INFO-all --stdout | bgzip -c > $FILTERED_VCF

# as good practice, we should index this too
bcftools index $FILTERED_VCF

# have a look at the log file
cat out.log

# next we will run plink
cd ../plink

# linkage pruning with plink
plink --vcf $FILTERED_VCF --double-id --allow-extra-chr --chr-set 30 \
--set-missing-var-ids @:# \
--maf 0.10 --indep-pairwise 50 10 0.3 --out sparrows_maf10

# extract the actual data
plink --vcf $FILTERED_VCF --double-id --allow-extra-chr --chr-set 30 \
--set-missing-var-ids @:# \
--extract sparrows_maf10.prune.in --make-bed --out sparrows_maf10_ld_prune

--out 