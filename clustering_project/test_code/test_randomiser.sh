#!/bin/bash

export TEMPORARY_DIRECTORY=`mktemp -d /tmp/XXXXXX`
python test_code/test_randomiser.py $1 $2 $3
rm -r $TEMPORARY_DIRECTORY
unset TEMPORARY_DIRECTORY
