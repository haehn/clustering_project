#!/bin/bash

echo 'Arguments: 1=directory suffix($1), 2=number of NNIs($2), 3=program($3), 4=model($4), 5=ncat($5), 6=datatype($6), 7=regime($7)'

WRITEDIR=/ebi/research/goldman/kevin/phymlreg3/sim$1
TMPWRAP=/net/isilon7/nobackup/research/goldman/kevin/clustering_project/tempdir_wrapper.sh
SIMULATE=/net/isilon7/nobackup/research/goldman/kevin/clustering_project/simulate.py
PICKLE=/net/isilon7/nobackup/research/goldman/kevin/clustering_test/pickle_simulation.py
ANALYSE=/net/isilon7/nobackup/research/goldman/kevin/clustering_test/analyse.py

mkdir $WRITEDIR
python $SIMULATE -k 4 -m 60 -n 16 -d $WRITEDIR -tmp $TEMPORARY_DIRECTORY -ratevar -nni $2 -r $7 -q -g $GTP_PATH
python $PICKLE -d $WRITEDIR -g $GTP_PATH -p $3 -m $4 -n $5 -data $6

bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -c ward     -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -c single   -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -c average  -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -c complete -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -c spectral -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -c MDS      -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d geo -c kmedoids -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -c ward     -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -c single   -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -c average  -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -c complete -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -c spectral -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -c MDS      -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d euc -c kmedoids -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -c ward     -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -c single   -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -c complete -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -c average  -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -c spectral -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -c MDS      -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
bsub -q research-rh6 -o /dev/null bash $TMPWRAP python $ANALYSE -d sym -c kmedoids -i $WRITEDIR -p $3 -m $4 -n $5 -data $6
