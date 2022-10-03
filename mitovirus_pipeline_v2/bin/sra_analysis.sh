#!/bin/bash

codonValues="4"
source_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P)

echo $source_dir
touch mitovirus_summary.csv

for codon in $codonValues
do
    mkdir "blastxcodontable"$codon
    cd "blastxcodontable"$codon
    diamond blastx -d $1 -q ../$2 --query-gencode $codon --sensitive --min-orf 300 -o ./blastx_output.tsv
    mkdir virusfams
    python $source_dir/AlignmentBreakup.py -i blastx_output.tsv
    mv *.txt virusfams
    # here we subset the nodes into their families
    $source_dir/node_split.sh ../$2 $codon
    # pull together all summaries
    for f in virusfams/*
    do
	if [[ $f == 'virusfams/Mitovirus' ]]
	then
	    cat virusfams/Mitovirus/Mitovirus_summary.csv >> ../mitovirus_summary.csv
	fi
    done   
    cd ..
done
