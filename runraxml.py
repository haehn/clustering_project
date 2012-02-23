#!/usr/bin/env python

import glob,os,shlex,subprocess,sys
from raxml_wrapper import run_RAxML
from handleArgs import handleArgs

command_line_args = handleArgs(sys.argv,help='''
runraxml.py arguments:
         -l = Location, either 'cluster' or 'local' (default) 
''')

besttree_path = "./trees/besttrees"
info_path = "./trees/info"
if not os.path.exists(besttree_path): os.mkdir(besttree_path)
if not os.path.exists(info_path): os.mkdir(info_path)

MSAs = glob.glob('./MSA/*.phy') #all per-locus *.phy alignment files in a list
for alignment in MSAs:
    name = 'gene'+alignment[alignment.index('_')+1:alignment.rindex('.')] #RAxML appends a name to the output of a job, which here is set to 'gene##'
    if command_line_args['-l']=="cluster": #define the shell command to run RAxML on the cluster, move the files we want to keep, and delete the rest
        process_args = 'bsub -o /dev/null -q research-rh6 "raxml -T 4 -m PROTGAMMAWAG -s {0} -n {1} -N 10 && mv RAxML_bestTree.{1} ./trees/besttrees/{1}.nwk && mv RAxML_info.{1} ./trees/info/{1}.info && rm *.{1} "'.format(alignment, name)
        subprocess.call( shlex.split(process_args) )
    elif command_line_args['-l']=="local": #define the shell command to run RAxML locally, move the files we want to keep, and delete the unwanted files
        process_args = 'raxml -T 4 -m PROTGAMMAWAG -s {0} -n {1} && mv RAxML_bestTree.{1} ./trees/besttrees/{1}.nwk && mv RAxML_info.{1} ./trees/info/{1}.info && rm *.{1} '.format(alignment,name)
        subprocess.call( shlex.split(process_args) )