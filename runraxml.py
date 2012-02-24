#!/usr/bin/env python

import glob,os,shlex,subprocess,sys
from handleArgs import handleArgs

command_line_args = handleArgs(sys.argv,help='''
runraxml.py arguments:
         -l = Location, either 'cluster' or 'local' (default)
         -dna = invokes GTRGAMMA model in RAxML (default is PROTGAMMAWAG)
         -dir = input directory 
''')

#Check command line arguments
if not '-dir' in command_line_args:
    INPUT_DIR = '.'
elif not command_line_args['-dir']:
    print "No input directory specified, trying '.'"
    INPUT_DIR = '.'
else: INPUT_DIR = command_line_args['-dir']

if '-dna' in command_line_args: 
    model = 'GTRGAMMA'
else: model = 'PROTGAMMAWAG'

besttree_path = "{}/trees/besttrees".format(INPUT_DIR)
info_path = "{}/trees/info".format(INPUT_DIR)
if not os.path.exists(besttree_path): os.mkdir(besttree_path)
if not os.path.exists(info_path): os.mkdir(info_path)

MSAs = glob.glob('{}/MSA/*.phy'.format(INPUT_DIR)) #all per-locus *.phy alignment files in a list
for alignment in MSAs:
    name = alignment[alignment.rindex('/')+1:alignment.rindex('.')] #RAxML appends a name to the output of a job, which here is set to 'gene##'
    #name = alignment[alignment.rindex('/')+1:alignment.rindex('.')]
    if command_line_args['-l']=="cluster": #define the shell command to run RAxML on the cluster, move the files we want to keep, and delete the rest
        process_args = 'bsub -o /dev/null -q research-rh6 "raxml -T 4 -m {0} -s {1} -n {2} -N 10 && mv RAxML_bestTree.{2} {3}/trees/besttrees/{2}.nwk && mv RAxML_info.{2} {3}/trees/info/{2}.info && rm *.{2} "'.format(model, alignment, name, INPUT_DIR)
        subprocess.call( shlex.split(process_args) )
    elif command_line_args['-l']=="local": #define the shell command to run RAxML locally, move the files we want to keep, and delete the unwanted files
        os.system( 'raxml -T 8 -m {0} -s {1} -n {2} -p 121 && mv RAxML_bestTree.{2} {3}/trees/besttrees/{2}.nwk && mv RAxML_info.{2} {3}/trees/info/{2}.info && rm *.{2} '.format(model, alignment,name,INPUT_DIR) )
        
#Why os.system works on the local call above and not subprocess.call - I DON'T KNOW
