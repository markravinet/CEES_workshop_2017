#!/bin/bash

## define variables
PICARD_PATH=/projects/cees/bin/picard/2.7.1/picard.jar
GATK_PATH=/cluster/software/VERSIONS/gatk-3.7/GenomeAnalysisTK.jar
REF=/work/users/msravine/ref/house_sparrow_genome_assembly-18-11-14_masked.fa
GVCF_LIST=/work/users/msravine/course_reads/sparrow_gvcf.list
OUTPUT_VCF=/work/users/msravine/course_reads/vcf/house_spanish_italian_toy.vcf

# run joint genotyping
echo "Running joint genotyping for gvcfs"
java -Xmx2g -jar $GATK_PATH -T GenotypeGVCFs \
-R $REF \
-V $GVCF_LIST \
-o $OUTPUT_VCF
