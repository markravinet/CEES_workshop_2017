#!/bin/sh

module load samtools bwa/0.7.10

INPUT_DIR=/work/users/msravine/workshop/trimmed
OUTPUT_DIR=/work/users/msravine/workshop/bams
REF=/projects/cees2/in_progress/cees_workshop/ref/house_sparrow_genome_assembly-18-11-14_masked.fa
SAMPLES=(house_8L19786 italian_Guglionesi_426 spanish_Lesina_285)

for SAMPLE in ${SAMPLES[@]}
do
F_PAIR=${INPUT_DIR}/${SAMPLE}_R1_trim_pair.fastq.gz
R_PAIR=${INPUT_DIR}/${SAMPLE}_R2_trim_pair.fastq.gz
F_UNPAIR=${INPUT_DIR}/${SAMPLE}_R1_trim_unpair.fastq.gz
R_UNPAIR=${INPUT_DIR}/${SAMPLE}_R2_trim_unpair.fastq.gz

# echo $F_PAIR $R_PAIR $F_UNPAIR $R_UNPAIR
# echo "Aligning $F_PAIR & $R_PAIR"
bwa mem $REF $F_PAIR $R_PAIR | samtools view -b | samtools sort > ${OUTPUT_DIR}/${SAMPLE}_pair.bam

bwa mem $REF $F_UNPAIR | samtools view -b | samtools sort > ${OUTPUT_DIR}/${SAMPLE}_f_unpair.bam

bwa mem $REF $R_UNPAIR | samtools view -b | samtools sort > ${OUTPUT_DIR}/${SAMPLE}_r_unpair.bam

samtools merge -f ${OUTPUT_DIR}/${SAMPLE}_merge.bam ${OUTPUT_DIR}/${SAMPLE}_pair.bam \
${OUTPUT_DIR}/${SAMPLE}_f_unpair.bam \
${OUTPUT_DIR}/${SAMPLE}_r_unpair.bam

samtools sort ${OUTPUT_DIR}/${SAMPLE}_merge.bam > ${OUTPUT_DIR}/${SAMPLE}_merge_sort.bam

done