#!/bin/bash

## define variables
PICARD_PATH=/projects/cees/bin/picard/2.7.1/picard.jar
GATK_PATH=/cluster/software/VERSIONS/gatk-3.7/GenomeAnalysisTK.jar
REF=/work/users/msravine/ref/house_sparrow_genome_assembly-18-11-14_masked.fa

# define directories
INPUT_DIR=/work/users/msravine/course_reads/bamfiles_gatk
OUTPUT_DIR=/work/users/msravine/course_reads/gvcf

# generate an array of samples
SAMPLES=($(for bam in ${INPUT_DIR}/*_realigned.bam; do echo $(basename ${bam%_*}); done))

for SAMPLE in ${SAMPLES[@]}
do

## RUN GATK HAPLOTYPE CALLER ###
# perform Haplotype calling
echo "### Calling variants for $SAMPLE ###"
java -Xmx2g -jar $GATK_PATH -T HaplotypeCaller \
-R $REF \
-I ${INPUT_DIR}/${SAMPLE}_realigned.bam \
--emitRefConfidence GVCF \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-o ${OUTPUT_DIR}/${SAMPLE}_raw.g.vcf

done