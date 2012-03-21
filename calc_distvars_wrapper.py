#!/usr/bin/env python

import os, sys

fil = sys.argv[1]
seqtype = sys.argv[2]
name = fil[:fil.rindex('.')]

os.system( "echo \"fil := ReadFastaWithNames('{0}'); seqtype := '{1}'; name := '{2}'; ReadProgram('calc_distvars.drw'); \" | darwin ".format(fil, seqtype, name) )