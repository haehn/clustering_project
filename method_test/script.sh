#!/bin/bash

# Makes a number of simulated data sets using LSB JobArrays
# Usage: bsub -o /dev/null (-e error.%J) -J "jobname[1-100]" bash ~/research/clustering_project/tempdir_wrapper.sh bash script.sh ARGS

echo 'Arguments: 1=number of classes, 2=number of species, 3=number of genes, 4=regime,'
echo '           5=master tree generator 6=permutation type 7=permutation level,'
echo number of classes  = $1
echo number of species  = $2
echo number of genes    = $3
echo regime             = $4
echo master tree gen    = $5
echo permutation type   = $6
echo permutation level  = $7

WRITEDIR=/ebi/research/goldman/kevin/regime$4/$6/level$7/sim$LSB_JOBINDEX
TMPWRAP=/nfs/nobackup/research/goldman/kevin/clustering_project/tempdir_wrapper.sh
SIMULATE=/nfs/nobackup/research/goldman/kevin/clustering_project/simulate.py

mkdir -p $WRITEDIR
python $SIMULATE -k $1 -n $2 -m $3 -r $4 -t 0 -master=$5 -c=$6 -p $7 -d $WRITEDIR -tmp $TEMPORARY_DIRECTORY -ratevar -q -u -g $GTP_PATH
