#!/bin/bash

WD1=/scratch/rtraborn/Daphnia/PA42_genome
GENOME=Daphnia_pulex_PA42_v3.0.fasta
WD2=/scratch/rtraborn/Daphnia/CAGE_reads/demultiplexed_fastqs
DEST_DIR=/scratch/rtraborn/Daphnia/CAGE_reads/aligned_bams

cd $WD1

echo "Indexing the genome"

bwa index $GENOME

cd $WD2

echo "Aligning the CAGE reads to PA42v3"

for FQ in *.fastq
do
bwa aln -t8 -B 3 -n 3 $GENOME -f $(basename $FQ .fastq).sai $FQ
bwa samse $GENOME $(basename $FQ .fastq).sai $FQ |
    samtools view -uS - |
    samtools sort - $(basename $FQ .fastq)
done

mv *.bam $DEST_DIR

echo "Job Complete!"

