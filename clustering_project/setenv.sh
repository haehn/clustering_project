#!/bin/bash

SCRIPT_PATH="${BASH_SOURCE[0]}";
if ([ -h "${SCRIPT_PATH}" ]) then
  while([ -h "${SCRIPT_PATH}" ]) do SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
pushd . > /dev/null
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`/bin/pwd`;
popd  > /dev/null

if ! $(echo "$PYTHONPATH" | tr ":" "\n" | grep -qx "$SCRIPT_PATH/class_files") ; then export PYTHONPATH=$SCRIPT_PATH/class_files:$PYTHONPATH ; fi
export DARWINHELPER=$SCRIPT_PATH/class_files/DV_wrapper.drw
