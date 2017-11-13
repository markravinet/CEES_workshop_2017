#!/bin/sh

# set up modules
module purge
module load bwa/0.7.10
module load samtools

## set variables ##
# base directory for files
OUTDIR=/work/users/msravine/course_reads/bamfiles

# set reference
REF=/work/users/msravine/ref/house_sparrow_genome_assembly-18-11-14_masked.fa 

cd /work/users/msravine/course_reads/full_trimmed

PREFIXES=($(for i in *_R1_trim_pair.fastq.gz; do echo ${i%%_R*}; done | uniq))

for PREFIX in ${PREFIXES[@]}
do

# paired
FORWARD=/work/users/msravine/course_reads/full_trimmed/${PREFIX}_R1_trim_pair.fastq.gz
REVERSE=/work/users/msravine/course_reads/full_trimmed/${PREFIX}_R2_trim_pair.fastq.gz
#unpaired
FORWARD_UNPAIR=/work/users/msravine/course_reads/full_trimmed/${PREFIX}_R1_trim_unpair.fastq.gz
REVERSE_UNPAIR=/work/users/msravine/course_reads/full_trimmed/${PREFIX}_R2_trim_unpair.fastq.gz

# nn9244k
# remove read read pair suffix
PREFIX=${PREFIX%_R*}
# get sample ID
SAMPLE=${PREFIX#*.*.*.*.*} 
# get read group ID
ID=${PREFIX%.*.*}
# get library
LIBRARY=${PREFIX%.*}
LIBRARY=${LIBRARY##*.}
# set platform
PLATFORM=Illumina
# extract PU data
PU_DATA=$(zcat $FORWARD | head -1 | cut -d ":" -f 3,4)

# construct read group
READGROUP="@RG\tID:${ID}\tPL:${PLATFORM}\tLB:${LIBRARY}\tSM:${SAMPLE}\tPU:${PU_DATA}"

# set ID
ID=$(basename $FORWARD)
ID=${ID%%_*}
echo $ID

# working directory
# set WORKDIR
WORKDIR=/work/users/msravine/course_reads/bamfiles

# move to workdir
cd $WORKDIR

# map paired end
echo "Aligning $SAMPLE, paired-end"
bwa mem -M -t 16 \
-R $READGROUP \
$REF $FORWARD $REVERSE | samtools view -b | samtools sort -T ${SAMPLE} > ${WORKDIR}/${SAMPLE}_pe.bam

# map unpaired - FORWARD
echo "Aligning $SAMPLE, unpaired, forward"
bwa mem -M -t 16 \
-R $READGROUP \
$REF $FORWARD_UNPAIR | samtools view -b | samtools sort -T ${SAMPLE} > ${WORKDIR}/${SAMPLE}_up_forward.bam

# map unpaired - REVERSE
echo "Aligning $SAMPLE, unpaired, reverse"
bwa mem -M -t 16 \
-R $READGROUP \
$REF $REVERSE_UNPAIR | samtools view -b | samtools sort -T ${SAMPLE} > ${WORKDIR}/${SAMPLE}_up_reverse.bam

## merge paired and unpaired
echo "Merging $SAMPLE bams"
samtools merge -rf ${WORKDIR}/${SAMPLE}_merge.bam \
${WORKDIR}/${SAMPLE}_pe.bam  ${WORKDIR}/${SAMPLE}_up_forward.bam ${WORKDIR}/${SAMPLE}_up_reverse.bam

# sort merged
echo "Sorting merged $SAMPLE bam"
samtools sort -T ${SAMPLE} -o ${WORKDIR}/${SAMPLE}_merge_sort.bam ${WORKDIR}/${SAMPLE}_merge.bam

# rm unsorted
rm ${WORKDIR}/${SAMPLE}_merge.bam ${WORKDIR}/${SAMPLE}_pe.bam ${WORKDIR}/${SAMPLE}_up_forward.bam ${WORKDIR}/${SAMPLE}_up_reverse.bam

done

