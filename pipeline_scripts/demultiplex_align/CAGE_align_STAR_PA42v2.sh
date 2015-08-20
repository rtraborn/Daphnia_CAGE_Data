#!/bin/bash

GENOME=/home/rtraborn/Daphnia/genome_sequence/PA42_v2_assembly/PA42_scaffold2.0.fasta
GENDIR=/home/rtraborn/Daphnia/genome_sequence/PA42_v2_assembly
WD=/home/rtraborn/Daphnia/CAGE/PA42_v2/demultiplexed_matched

cd $WD

#echo "Indexing the PA42 genome"

#STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $GENDIR --genomeFastaFiles $GENOME

echo "Aligning the CAGE reads to PA42v2"

echo "Mature females"

#cd $WD/mat_females/mat_females_1

#    STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_mat_females_1.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10

#cd $WD/mat_females/mat_females_2

#STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_mat_females_2.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10

cd $WD/mat_females/mat_females_3

    STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_mat_females_3.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10

echo "Mature Males"
    
cd $WD/mat_males/mat_males_1
    
    STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_mature_males_1.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10

cd $WD/mat_males/mat_males_2
    
    STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_mature_males_2.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10

echo "pE Females"
    
cd $WD/pE_females/pE_females_1
    
STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_pE_females_1.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10

cd $WD/pE_females/pE_females_2

STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_pE_females_2.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10

cd $WD/pE_females/pE_females_3

    STAR --runThreadN 8 --runMode alignReads --genomeDir $GENDIR -sjdbOverhang 46 --readFilesIn $WD/Daphnia_pE_females_3.fastq --readMapNumber -1 --clip5pNbases 3 --outStd BAM_SortedByCoordinate --outSAMtype BAM Unsorted --outSAMprimaryFlag AllBestScore --outFilterType Normal --outFilterMismatchNoverLmax 0.10
    
echo "Job Complete!"

