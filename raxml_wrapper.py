#!/usr/bin/env python

def run_RAxML(executable, sequence_file, model, run_name, nruns=1, random_tree=False, pthreads = False, bsub=False, nthreads = 2):
    import subprocess, shlex
    if bsub: executable = 'bsub -o /dev/null {0}'.format(executable)
    if pthreads: args = "{0} -T {1} -s {2} -m {3} -n {4} -N {5}".format(executable, str(nthreads), sequence_file, model, run_name, str(nruns))
    else: args = args = "{0} -s {1} -m {2} -n {3} -N {4}".format(executable, sequence_file, model, run_name, str(nruns))
    if random_tree: args = args + ' -d'
    subprocess.call(shlex.split(args))
