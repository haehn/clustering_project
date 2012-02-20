#!/usr/bin/env python

import glob,os

def run_RAxML(executable, sequence_file, model, run_name, nruns=1, random_tree=False, pthreads = False, bsub=False, nthreads = 2):
    import subprocess
    if bsub: executable = 'bsub -o /dev/null {0}'.format(executable)
    if pthreads: args = (executable, '-T', str(nthreads), '-s', sequence_file, '-m', model, '-n', run_name, '-N', str(nruns))
    else: args = (executable, '-s', sequence_file, '-m', model, '-n', run_name, '-N', nruns)
    if random_tree: args = args + ('-d',)
    subprocess.call(args)

besttree_path = "./trees/besttrees"
info_path = "./trees/info"
if not os.path.exists(besttree_path): os.mkdir(besttree_path)
if not os.path.exists(info_path): os.mkdir(info_path)

MSAs = glob.glob('./MSA/*.phy')
for alignment in MSAs[0:1]:
    name = 'gene'+alignment[alignment.index('_')+1:alignment.rindex('.')]
    run_RAxML('raxml', alignment, 'PROTGAMMAWAG', name, nruns=10, 
    random_tree=True, pthreads=True, nthreads=8, bsub=True)
    os.rename('./RAxML_bestTree.{0}'.format(name), './trees/besttrees/bestTree.{0}'.format(name))
    os.rename('./RAxML_info.{0}'.format(name), './trees/info/info.{0}'.format(name))
