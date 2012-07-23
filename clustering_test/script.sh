#!/bin/bash

WRITEDIR=/ebi/research/goldman/kevin/sim$1
TMPWRAP=/net/isilon7/nobackup/research/goldman/kevin/clustering_project/tempdir_wrapper.sh
SIMULATE=/net/isilon7/nobackup/research/goldman/kevin/clustering_project/simulate.py
PICKLE=/net/isilon7/nobackup/research/goldman/kevin/clustering_test/pickle_simulation.py
ANALYSE=/net/isilon7/nobackup/research/goldman/kevin/clustering_test/analyse.py


mkdir $WRITEDIR
#python $SIMULATE -k 4 -m 60 -n 16 -d $WRITEDIR -tmp $TEMPORARY_DIRECTORY -ratevar -nni $2 -r 2 -q
python $PICKLE -d $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -m ward -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -m single -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -m average -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -m complete -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -m spectral -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -m MDS -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -m kmedoids -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -m ward -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -m single -i $WRITEDIR
bsub -q research-rh6 bash $TMPWRAP python $ANALYSE -d euc -m average -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -m complete -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -m spectral -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -m MDS -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -m kmedoids -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -m ward -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -m single -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -m complete -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -m average -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -m spectral -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -m MDS -i $WRITEDIR
#bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -m kmedoids -i $WRITEDIR
