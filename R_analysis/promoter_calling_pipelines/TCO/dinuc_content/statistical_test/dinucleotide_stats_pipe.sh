for i in {1..100}; do

    bedtools slop -l $i -r -$i -s -i $1 -g ~/Daphnia/annotation_files/Scaffold_size.txt > dinuc_${i}.bed
    bedtools getfasta -s -fi /home/rtraborn/Daphnia/genome_sequence/TCO_assembly/Daphnia_pulex.fasta -bed dinuc_${i}.bed -fo dinuc_${i}.fasta
    homerTools freq -format fasta dinuc_${i}.fasta > dinuc_frequency_tables/dinuc_${i}.freq
    rm dinuc_${i}.bed dinuc_${i}.fasta
done
