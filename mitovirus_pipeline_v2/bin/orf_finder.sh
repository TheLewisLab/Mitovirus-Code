#!/bin/bash

source_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P)

for dir in */
do
    cd $dir
    #rm node_orfs.fa
    # find open reading frames using alternative start for codon table 2 and 5
    if [[ $1 == '4' || $1 == '1' || $1 == '3' ]]
    then
	$source_dir/ORFfinder -in *.fa -g $1 -s 0 -ml 400 -outfmt 1 -out node_orfs.fa
    else
	$source_dir/ORFfinder -in *.fa -g $1 -s 1 -ml 400 -outfmt 1 -out node_orfs.fa
    fi
    
    # run python codon frequency script
    python $source_dir/codonfrequencyanalysis.py -i node_orfs.fa -t $1 -o cds_codon_frequency.csv
    cd ..
done
