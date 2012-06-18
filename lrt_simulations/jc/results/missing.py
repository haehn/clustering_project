#!/usr/bin/env python

import glob,sys,pprint
import os
import time

if len(sys.argv)>3:
    if sys.argv[3]=='-launch':
        launch_jobs = True
    else:
        launch_jobs = False
else:
    launch_jobs = False

l = []
for i in range(1,20):
    for j in range(1,21):
        l.append( '{0}class_50sims_{2}_{1}'.format(i,j,sys.argv[2]))

files = glob.glob('{0}/*'.format(sys.argv[1]))

files = [x[x.rindex('/')+1:] for x in files if not 'classes' in x]

for f in files:
    if f in l:
        l.remove(f)
    else:
        print '{0} not in list'.format(f)
        sys.exit(0)

print l
if launch_jobs:
    for e in l:
        numclass = e.split('class')[0]
        numsim = e.split('_')[-1]
        os.system('bsub -q research-rh6 -M 8192 -R "rusage[mem=8192]" -o /dev/null bash tempdir_wrapper.sh python gc_script.py -o={0}/{1}class_50sims_{2}_{3} -b={4}_jc69_{5}.pickle -c={1} -s=50 -d={5}'.format(sys.argv[1],numclass,sys.argv[2],numsim,sys.argv[1].split('/')[-2],sys.argv[1].split('/')[-1]))
        time.sleep(0.5)
