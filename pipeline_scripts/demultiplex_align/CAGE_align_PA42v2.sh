#!/bin/bash

WD1=/home/rtraborn/Daphnia/genome_sequence/
GENOME=Daphnia_pulex.fasta
WD2=/home/rtraborn/Daphnia/CAGE/demultiplexed_matched/fastq_files

cd $WD1

echo "Indexing the genome"

bwa index $GENOME

cd $WD2

echo "Aligning the CAGE reads to TCO"

for FQ in *.fastq
do
    STAR --runThreadN 8 -runMode genomeGenerate --genomeDir /home/rtraborn/Daphnia/genome_annotation --genomeFastaFiles $GENOME --sjdbGTFfile $annotations -sjdbOverhang 46  -f $(basename $FQ .fastq).sai $FQ
bwa samse $GENOME $(basename $FQ .fastq).sai $FQ |
    samtools view -uS - |
    samtools sort - $(basename $FQ .fastq)
done

echo "Job Complete!"

