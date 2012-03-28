#!/bin/bash
DIR=real_data/yeast_data
#./runsim.py -c 4 -s 8 -dir $DIR
for i in {40..1}; do
    echo "./runquiet.py -dir $DIR/MSA -program treecollection -trees $DIR/chtr -clusters $DIR/chcl -m sym -l ward -nclasses $i >> $DIR/resultsrf.txt"
    ./runquiet.py -dir $DIR/MSA -program treecollection -trees $DIR/chtr -clusters $DIR/chcl -m sym -l ward -nclasses $i >> $DIR/results.txt
    echo "rm $DIR/chtr/cluster*"
    rm $DIR/chtr/cluster*
    echo "rm $DIR/chcl/*"
    rm $DIR/chcl/*
done
