#!/bin/bash

GENOME=/home/rtraborn/Daphnia/genome_sequence/Daphnia_pulex.fasta
WD=/home/rtraborn/Daphnia/CAGE/demultiplexed_matched/fastq_files
rRNA=/home/rtraborn/Daphnia/annotation_files/Daphnia_rDNA.fasta

cd $WD

echo "Aligning the CAGE reads to TCO"

for FQ in *.fastq
do
bwa aln -t16 -B 3 -n 3 $GENOME -f $(basename $FQ .fastq).sai $FQ
bwa samse $GENOME $(basename $FQ .fastq).sai $FQ |
    samtools view -uS - |
    samtools sort - $(basename $FQ .fastq) && samtools index $(basename $FQ .fastq).bam
done

echo "Removing rRNA contamination"

for BM in *.bam
do
echo $BM rRNAdust
  (rRNAdust -t16 $rRNA $BM -e 3 | samtools view -bS - 2> /dev/null | sponge $(basename $BM .bam).filtered.bam) 2>&1 | sed 's/Excluded: //'
  done | tee complete_rRNA_stats.log

echo "Job Complete!"
