#!/bin/bash

export TEMPORARY_DIRECTORY=`mktemp -d /tmp/XXXXXX`

python write_simulation_parameters.py $1

rm -r $TEMPORARY_DIRECTORY
unset TEMPORARY_DIRECTORY
