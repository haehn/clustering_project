#!/bin/bash

export TEMPORARY_DIRECTORY=`mktemp -d /tmp/XXXXXX`
python $1 $2 $3
rm -r $TEMPORARY_DIRECTORY
unset TEMPORARY_DIRECTORY
