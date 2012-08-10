#!/bin/bash

# Makes a number of simulated data sets using LSB JobArrays

echo 'Arguments: 1=directory suffix, 2=number of classes, 3=number of genes, 4=number of species, 5=number of NNIs,'
echo '           6=regime'
echo directory suffix  = $1
echo number of classes = $2
echo number of genes   = $3
echo number of species = $4
echo number of NNIs    = $5
echo regime            = $6

WRITEDIR=/ebi/research/goldman/kevin/regime2/nni$5/sim$LSB_JOBINDEX
TMPWRAP=/net/isilon7/nobackup/research/goldman/kevin/clustering_project/tempdir_wrapper.sh
SIMULATE=/net/isilon7/nobackup/research/goldman/kevin/clustering_project/simulate.py
PICKLE=/net/isilon7/nobackup/research/goldman/kevin/method_test/pickle_simulation.py
ANALYSE=/net/isilon7/nobackup/research/goldman/kevin/method_test/analyse.py

mkdir $WRITEDIR
python $SIMULATE -k $2 -m $3 -n $4 -d $WRITEDIR -tmp $TEMPORARY_DIRECTORY -ratevar -nni $5 --regime=$6 -q -g $GTP_PATH