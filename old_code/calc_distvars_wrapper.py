#!/usr/bin/env python

import os, sys

fil = sys.argv[1] # filename of fasta file to convert to DistVar matrix (.dv) file
seqtype = sys.argv[2] # either 'AA' or 'DNA'
calc_dv_drw_filepath = sys.argv[3] # path to calc_distvars.drw
name = fil[:fil.rindex('.')]

os.system( "echo \"fil := ReadFastaWithNames('{0}'); seqtype := '{1}'; name := '{2}'; ReadProgram('{3}'); \" | darwin ".format(fil, seqtype, name, calc_dv_drw_filepath) )
