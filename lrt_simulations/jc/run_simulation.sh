#!/bin/bash

export TEMPORARY_DIRECTORY=`mktemp -d /tmp/XXXXXX`

python run_simulation.py $1 $2

rm -r $TEMPORARY_DIRECTORY
unset TEMPORARY_DIRECTORY
