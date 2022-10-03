#!/bin/bash

touch Mitovirus_summary.csv

for node in ./*/
do
    cd $node
    cat summary.csv >> ../Mitovirus_summary.csv
    cd ..
done
