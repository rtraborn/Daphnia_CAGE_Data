#!/bin/bash

rRNA=/home/rtraborn/Daphnia/annotation_files/Daphnia_rDNA.fasta
for BM in *.bam
do
echo $BM rRNAdust
  (rRNAdust -t8 $rRNA $BM -e 3 | samtools view -bS - 2> /dev/null | sponge $(basename $BM .bam).rRNA.bam) 2>&1 | sed 's/Excluded: //'
  done | tee complete_rRNA_stats.log


