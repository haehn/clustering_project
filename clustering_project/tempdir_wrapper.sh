#!/bin/bash

export TEMPORARY_DIRECTORY=`mktemp -d /tmp/XXXXXXXX`

$*

rm -r $TEMPORARY_DIRECTORY
unset TEMPORARY_DIRECTORY
