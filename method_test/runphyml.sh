#!/bin/bash

# Run phyml on all *.phy files in a directory using LSB jobarrays
# Analysis options are hard-coded - why not?

echo 'Arguments: 1=Directory'

for F in $1/*.phy; do bsub -o output.%J -e error.%J phyml -i $F -m GTR -c 4 -d nt -a e -b 0 --sequential --no_memory_check --run_id ml; done
