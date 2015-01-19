GENOME=/home/rtraborn/Daphnia/genome_sequence/DP.fasta
for FQ in *.fastq
do
bwa aln -t6 -B 3 $GENOME -f $(basename $FQ .fastq).sai $FQ
bwa samse $GENOME $(basename $FS .fastq).sai $FQ
    samtools view -uS - |
    samtools sort - $(basename $FQ .fastq)
done
