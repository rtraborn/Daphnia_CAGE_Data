#!/bin/bash

WD1=/home/rtraborn/Daphnia/genome_sequence/PA42_v2_assembly/
GENOME=PA42_scaffold2.0.fasta 
WD2=/home/rtraborn/Daphnia/CAGE/demultiplexed_matched/fastq_files

cd $WD1

echo "Indexing the genome"

bwa index $GENOME

cd $WD2

echo "Aligning the CAGE reads to TCO"

for FQ in *.fastq
do
bwa aln -t12 -B 3 -n 3 $GENOME -f $(basename $FQ .fastq).sai $FQ
bwa samse $GENOME $(basename $FQ .fastq).sai $FQ |
    samtools view -uS - |
    samtools sort - $(basename $FQ .fastq)
done

echo "Job Complete!"

