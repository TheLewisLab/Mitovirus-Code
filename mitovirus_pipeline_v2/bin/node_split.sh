#!/bin/bash

source_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P)
bucketpath=$(pwd)

for f in virusfams/*
do
    echo $f
    # pull sequence fasta files for each node for each family
    seqtk subseq $1 $f > ${f%.*}.fa
    mkdir ${f%.*}/
    mv ${f%.*}.* ${f%.*}/
    cd ${f%.*}/
    # split by node
    while read line
    do
	if [[ ${line:0:1} == '>' ]]
	then
	    temp_arr=($line)
	    mkdir ${temp_arr[0]#>}/
	    outfile=${temp_arr[0]#>}/${temp_arr[0]#>}.fa
	    echo $line > $outfile
	else
	    echo $line >> $outfile
	fi
    done < *.fa
    # run orf finder module
    if [[ $f == 'virusfams/Mitovirus.txt' ]]
    then
	echo "inside mito"
	$source_dir/orf_finder.sh $2
	$source_dir/summarizer.sh
    fi
    cd $bucketpath
done
