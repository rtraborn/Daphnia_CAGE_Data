#!/bin/bash

#PBS -N Daphnia_rRNA_align_pipeline_mat_female1
#PBS -k o
#PBS -q cpu
#PBS -l nodes=1:ppn=32,vmem=200gb
#PBS -l walltime=8:00:00
#PBS -q shared
#BS -m abe
#PBS -M rtraborn@indiana.edu

WD=/N/dc2/projects/Daphnia_Gene_Expression/demultiplex/mature_females

ulimit -s unlimited
module load gcc/4.7.2
module load bwa/0.7.2
module load samtools
module load bedtools/2.20.1

cd $WD #change to the working directory

bwa mem -a -t 32 /N/dc2/projects/Daphnia_Gene_Expression/rRNA_seqs/Daphnia_rRNA.fasta demulti-mat_females_1.fastq > mat_females_1_rRNA.sam

samtools view -Sb mat_females_1_rRNA.sam > mat_females_1_rRNA.bam

bamToBed -i cigar mat_females_1_rRNA.bam > mat_females_1_rRNA.bed

