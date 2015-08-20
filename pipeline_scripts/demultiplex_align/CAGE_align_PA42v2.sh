#!/bin/bash

WD1=/home/rtraborn/Daphnia/genome_sequence
GENOME=Daphnia_pulex.fasta
WD2=/home/rtraborn/Daphnia/CAGE/PA42_v2/demultiplexed_matched

cd $WD1

echo "Indexing the genome"

bwa index $GENOME

cd $WD2

echo "Aligning the CAGE reads to PA42v2"

for FQ in *.fastq
do
    STAR --runThreadN 8 -runMode genomeGenerate --genomeDir $WD1 --genomeFastaFiles $GENOME --sjdbGTFfile - -sjdbOverhang 46 --readFilesIn $WD2 --readFilesCommand zcat - --readMapNumber -1 --clip5pNbases 3 --outFileNamePrefix ./ --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --outSAMprimaryFlag AllBestScore --outWigType bedGraph read1_5p --outWigStrand Stranded --outWigNorm RPM --outFilterType Normal --outFilterMismatchNoverLmax 0.2 
done

echo "Job Complete!"

