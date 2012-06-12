#!/usr/bin/env python

import re
import glob
import sys

directory = sys.argv[1]
files = glob.glob('{0}/*.pickle'.format(directory))
chk = range(1,101)

for f in files:
    chk.remove(int(re.search('(?<=result)(\d+)',f).group()))

print chk
 
