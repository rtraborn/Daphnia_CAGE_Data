for i in {1..10}; do
    findMotifsGenome.pl /home/rtraborn/Daphnia/CAGE/TCO/motif_analysis/Dp_homer_tag_directory/training_data/train_list_$i.txt ~/Daphnia/annotation_files/Daphnia_pulex.fasta homer_train_$i  -size -50,50 -len 6,8,10,12 -bits -mset insects
    done
    
