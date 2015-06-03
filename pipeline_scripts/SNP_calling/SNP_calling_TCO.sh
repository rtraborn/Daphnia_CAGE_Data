#!/bin/bash

WD=/home/rtraborn/Daphnia/CAGE/demultiplexed_matched/bam_files/TCO_bams/bam_filtered/files/

cd $WD

echo "Merging the alignment replicate files, as appropriate"

echo "Merging mature females"

samtools merge -r -@ 4 Dp_mat_females_filtered.bam Daphnia_mat_females_1.rRNA.bam Daphnia_mat_females_2.rRNA.bam Daphnia_mat_females_3.rRNA.bam

#echo "Merging pE females"

samtools merge -r -@ 4 Dp_pE_females_filtered.bam Daphnia_pE_females_1.rRNA.bam Daphnia_pE_females_2.rRNA.bam Daphnia_pE_females_3.rRNA.bam

echo "Merging mature males"

samtools merge -r -@ 4 Dp_mat_males_filtered.bam Daphnia_mature_males_1.rRNA.bam Daphnia_mature_males_2.rRNA.bam

echo "Calling SNPs from each merged alignment"

for BM in *_filtered.bam
do
  samtools mpileup -I -v -d 1000 -C 50 --output $(basename $BM .filtered.bam).vcf $BM
 
  done

echo "Job Complete!"

