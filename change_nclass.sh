#!/bin/bash
DIR=scratch/chclass$1
./runsim.py -c 4 -s 8 -dir $DIR
for i in {40..1}; do
    echo "runquiet.py -dir $DIR/MSA -program treecollection -trees $DIR/tr -clusters $DIR/cl -m sym -l ward -nclasses $i >> $DIR/results.txt"
    ./runquiet.py -dir $DIR/MSA -program treecollection -trees $DIR/tr -clusters $DIR/cl -m sym -l ward -nclasses $i >> $DIR/results.txt
    echo "rm $DIR/tr/cluster*"
    rm $DIR/tr/cluster*
    echo "rm $DIR/cl/*"
    rm $DIR/cl/*
done
