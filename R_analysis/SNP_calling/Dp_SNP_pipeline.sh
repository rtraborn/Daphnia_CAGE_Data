#/usr/bin/sh

#need to redo alingment to fix headers

GENOME=/home/rtraborn/Daphnia/genome_sequence/DP.fasta

samtools mpileup -uD -f $GENOME /home/rtraborn/Daphnia/CAGE/merged_bam/pE_fem_filtered_merged.bam /home/rtraborn/Daphnia/CAGE/merged_bam/mat_fem_filtered_merged.bam |
    bcftools view -bvcg - > Dp_filtered_merged.raw.bcf
