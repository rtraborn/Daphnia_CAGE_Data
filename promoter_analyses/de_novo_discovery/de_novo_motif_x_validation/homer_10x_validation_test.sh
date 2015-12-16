for i in {1..10}; do
    findMotifsGenome.pl /home/rtraborn/Daphnia/CAGE/TCO/motif_analysis/Dp_homer_tag_directory/test_data/test_list_$i.txt ~/Daphnia/annotation_files/Daphnia_pulex.fasta homer_test_$i -size -50,50 -len 6,8,10,12 -bits -mset insects
    done
    
