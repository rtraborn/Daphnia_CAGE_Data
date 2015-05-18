#!/bin/bash

FASTQ=Daphnia_CAGE_combined.fastq.gz

echo "Generatating Quality Report and Graph"

zcat $FASTQ | fastx_quality_stats -o $(basename $FASTQ .fastq.gz).fastx.qual
fastx_nucleotide_distribution_graph.sh -i $(basename $FASTQ .fastq.gz).fastx.qual -t "Daphnia CAGE Quality Stats" -o Daphnia_CAGE_Quality.png
