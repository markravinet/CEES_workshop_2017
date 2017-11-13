#!/bin/sh

for FASTQ in *.fastq.gz
do
echo $FASTQ
zcat $FASTQ | wc -l | awk '{print $1/4}'
done