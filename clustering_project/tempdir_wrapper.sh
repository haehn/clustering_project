#!/bin/bash

export TEMPORARY_DIRECTORY=`mktemp -d /tmp/tdwrap.XXXXXXXXXX`

$*

rm -r $TEMPORARY_DIRECTORY
unset TEMPORARY_DIRECTORY
