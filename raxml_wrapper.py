#!/usr/bin/env python

def run_RAxML(executable, sequence_file, model, run_name, pthreads = False, bsub=False, nthreads = 2):
    import subprocess
    if bsub: executable = 'bsub -o /dev/null {0}'.format(executable)
    if pthreads: args = (executable, '-T', str(nthreads), '-s', sequence_file, '-m', model, '-n', run_name)
    else: args = (executable, '-s', sequence_file, '-m', model, '-n', run_name)
    print args
    subprocess.call(args)

