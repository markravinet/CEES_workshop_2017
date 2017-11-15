#!/bin/bash

## define variables
PICARD_PATH=/projects/cees/bin/picard/2.7.1/picard.jar
GATK_PATH=/cluster/software/VERSIONS/gatk-3.7/GenomeAnalysisTK.jar
REF=/work/users/msravine/ref/house_sparrow_genome_assembly-18-11-14_masked.fa

# define directories
INPUT_DIR=/work/users/msravine/course_reads/bamfiles
OUTPUT_DIR=/work/users/msravine/course_reads/bamfiles_gatk

# generate an array of samples
SAMPLES=($(for bam in ${INPUT_DIR}/*.bam; do echo $(basename ${bam%.*}); done))

#SAMPLE=${SAMPLES[0]}

for SAMPLE in ${SAMPLES[@]}
do

# first sort using PICARD
echo "### Running Picard sort on $SAMPLE ###"
java -Xmx2g -jar $PICARD_PATH SortSam \
INPUT=$INPUT_DIR/${SAMPLE}.bam \
OUTPUT=${OUTPUT_DIR}/${SAMPLE}_sort.bam \
SORT_ORDER=coordinate

# then mark duplicates
echo "### Running Picard MarkDuplicates on $SAMPLE ###"
java -Xmx2g -jar $PICARD_PATH MarkDuplicates \
INPUT=${OUTPUT_DIR}/${SAMPLE}_sort.bam \
OUTPUT=${OUTPUT_DIR}/${SAMPLE}_sort_dedup.bam \
METRICS_FILE=${OUTPUT_DIR}/${SAMPLE}_dedup_metrics.txt 

# index bams using PICARD
echo "### Running Picard sort on $SAMPLE ###"
java -Xmx2g -jar $PICARD_PATH BuildBamIndex \
INPUT=${OUTPUT_DIR}/${SAMPLE}_sort_dedup.bam

## RUN GATK MODULES ###
# find indel realignment intervals
echo "### Identifying realignment targets for $SAMPLE ###"
java -Xmx2g -jar $GATK_PATH -T RealignerTargetCreator \
-R $REF \
-I ${OUTPUT_DIR}/${SAMPLE}_sort_dedup.bam \
-o ${OUTPUT_DIR}/${SAMPLE}_sort_dedup.intervals \
-nt 2

# perform indel realignment
echo "### Performing realignment for $SAMPLE ###"
java -Xmx2g -jar $GATK_PATH -T IndelRealigner \
-R $REF \
-I ${OUTPUT_DIR}/${SAMPLE}_sort_dedup.bam \
-targetIntervals ${OUTPUT_DIR}/${SAMPLE}_sort_dedup.intervals \
-o ${OUTPUT_DIR}/${SAMPLE}_realigned.bam

# use Picard to index realigned bam
echo "### Running Picard index on $SAMPLE ###"
java -Xmx2g -jar $PICARD_PATH BuildBamIndex \
INPUT=${OUTPUT_DIR}/${SAMPLE}_realigned.bam

done