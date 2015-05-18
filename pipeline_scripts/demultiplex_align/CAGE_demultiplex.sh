#!/bin/bash

BARCODES=barcodes.txt
FASTQ=Daphnia_CAGE_combined.fastq.gz

zcat $FASTQ | fastx_barcode_splitter.pl --bcfile $BARCODES --prefix Daphnia_ --suffix .fastq --bol --exact | sed 1d | cut -f1,2 | perl -ne 'print "extracted\t$_"' | grep --color=auto -v -e unmatched -e total | tee Daphnia_CAGE.extracted.log 
