#!/usr/bin/env python

import glob,os,sys
from raxml_wrapper import run_RAxML
from handleArgs import handleArgs

args = handleArgs(sys.argv,help='''
runraxml.py arguments:
         -l = Location, either 'cluster' or 'local' (default) 
''')

besttree_path = "./trees/besttrees"
info_path = "./trees/info"
if not os.path.exists(besttree_path): os.mkdir(besttree_path)
if not os.path.exists(info_path): os.mkdir(info_path)

MSAs = glob.glob('./MSA/*.phy')
for alignment in MSAs:
    name = 'gene'+alignment[alignment.index('_')+1:alignment.rindex('.')]
    print alignment, name
    if args['-l']=="cluster":
        run_RAxML('raxml', alignment, 'PROTGAMMAWAG', name, nruns=10, 
    random_tree=True, pthreads=True, nthreads=8, bsub=True)
    elif args['-l']=="local":
        run_RAxML('raxml', alignment, 'PROTGAMMAWAG', name, nruns=1, 
    random_tree=False, pthreads=False, bsub=False)