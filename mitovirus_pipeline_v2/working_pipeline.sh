#!/bin/bash

source_dir=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P)

echo $source_dir

for dir in $1*
do
    cd $dir
    # run sra module with $2 being protein ref
    $source_dir/bin/sra_analysis.sh $2 *.fasta
    cd ..
done
