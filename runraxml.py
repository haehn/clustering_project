#!/usr/bin/env python

import glob,os,shlex,subprocess,sys
from raxml_wrapper import run_RAxML
from handleArgs import handleArgs

command_line_args = handleArgs(sys.argv,help='''
runraxml.py arguments:
         -l = Location, either 'cluster' or 'local' (default) 
''')

MSAs = glob.glob('./MSA/*.phy')
for alignment in MSAs[0:1]:
    name = 'gene'+alignment[alignment.index('_')+1:alignment.rindex('.')]
    print alignment, name
    if command_line_args['-l']=="cluster":
        process_args = 'bsub -o /dev/null -q research-rh6 "raxml -T 4 -m PROTGAMMAWAG -s {0} -n {1} -N 10 && mv RAxML_bestTree.{1} ./trees/besttrees/{1}.nwk && mv RAxML_info.{1} ./trees/info/{1}.info && rm RAxML* "'.format(alignment, name)
        subprocess.call( shlex.split(process_args) )
    elif command_line_args['-l']=="local":
        process_args = 'raxml -T 4 -m PROTGAMMAWAG -s {0} -n {1} && mv RAxML_bestTree.{1} ./trees/besttrees/{1}.nwk && mv RAxML_info.{1} ./trees/info/{1}.info && rm RAxML* '.format(alignment,name)
        subprocess.call( shlex.split(process_args) )