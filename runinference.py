#!/usr/bin/env python

import glob,os,shlex,subprocess,sys
from handleArgs import handleArgs

def run_raxml(model, alignment, name, working_directory, tree_directory, info_directory): 
    os.system( 'raxml -T 4 -m {0} -s {1} -n {2} -p 121 > /dev/null && mv RAxML_bestTree.{2} {3}/{2}.nwk > /dev/null && mv RAxML_info.{2} {4}/{2}.info && rm *.{2} '.format(model, alignment, name, tree_directory, info_directory) )

def run_phyml(model, alignment, name, working_directory, tree_directory, info_directory, datatype):
    os.system( 'phyml -m {0} -i {1} -d {6} -a e -b 0 --sequential > /dev/null && mv {1}_phyml_tree.txt {4}/{2}.nwk && mv {1}_phyml_stats.txt {5}/{2}.info'.format( model, alignment, name, working_directory, tree_directory, info_directory, datatype))

def run_treecollection(alignment, name, tree_directory, info_directory, datatype):
    os.system( 'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\' ; ReadProgram(\'TC_wrapper.drw\');" | darwin > /dev/null'.format(alignment, datatype))
    os.system( 'TreeCollection temp_distvar.txt temp_map.txt temp_labels.txt temp_tree.nwk {0}/{1}.nwk > {2}/{1}.txt'.format(tree_directory, name, info_directory))
    os.system( 'rm temp_*')

##################################################################################################

command_line_args = handleArgs(sys.argv,help='''
runml.py arguments:
         -program = raxml or phyml
         -dna     = invokes GTRGAMMA model in RAxML (default is PROTGAMMAWAG)
         -dir     = input directory 
''')

#Check command line arguments

if not '-program' in command_line_args:
    print "No program chosen"
    sys.exit(1)
elif not command_line_args['-program']:
    print "No program chosen"
    sys.exit(1)
else:
    if command_line_args['-program'] == 'raxml':
        PROG_INDEX = 0
    elif command_line_args['-program'] == 'phyml':
        PROG_INDEX = 1
    elif command_line_args['-program'] == 'treecollection':
        PROG_INDEX = 2
    else:
        print "Unrecognised program chosen - use 'raxml', 'phyml' or 'treecollection'"
        sys.exit(1)

if not '-dir' in command_line_args:
    INPUT_DIR = '.'
elif not command_line_args['-dir']:
    print "No input directory specified, trying '.'"
    INPUT_DIR = '.'
else: INPUT_DIR = command_line_args['-dir']

if '-dna' in command_line_args: 
    raxml_model = 'GTRGAMMA'
    phyml_model = 'GTR'
    datatype = 'nt'
else: 
    raxml_model = 'PROTGAMMAWAG'
    phyml_model = 'WAG'
    datatype = 'aa'

##################################################################################################

treepaths = ["{0}/trees/raxmltrees".format(INPUT_DIR), "{0}/trees/phymltrees".format(INPUT_DIR), "{0}/trees/TCtrees".format(INPUT_DIR)]
infopaths = ["{0}/trees/raxmlinfo".format(INPUT_DIR), "{0}/trees/phymlinfo".format(INPUT_DIR), "{0}/trees/TCinfo".format(INPUT_DIR)] # No info yet for TC

if not os.path.exists("{0}/trees".format(INPUT_DIR)): os.mkdir("{0}/trees".format(INPUT_DIR))
if not os.path.exists(treepaths[PROG_INDEX]): os.mkdir(treepaths[PROG_INDEX])
if not os.path.exists(infopaths[PROG_INDEX]): os.mkdir(infopaths[PROG_INDEX])

MSAs = glob.glob('{}/MSA/*.phy'.format(INPUT_DIR)) #all per-locus *.phy alignment files in a list
models = [raxml_model, phyml_model, None]
for alignment in MSAs:
    model = models[PROG_INDEX]
    name = alignment[alignment.rindex('/')+1:alignment.rindex('.')] #RAxML appends a name to the output of a job, which here is set to 'gene##'
    tree_directory = treepaths[PROG_INDEX]
    info_directory = infopaths[PROG_INDEX]
    if PROG_INDEX == 0:
        print "Running raxml on {0}".format(name)
        run_raxml(model, alignment, name, INPUT_DIR, tree_directory, info_directory)
    elif PROG_INDEX == 1:
        print "Running phyml on {0}".format(name)
        run_phyml(model, alignment, name, INPUT_DIR, tree_directory, info_directory, datatype)
    elif PROG_INDEX == 2:
        if datatype == 'aa': datatype = 'AA'
        elif datatype == 'nt': datatype = 'DNA'
        alignment = alignment.replace(".phy", ".fas")
        print "Running TreeCollection on {0}".format(name)
        run_treecollection(alignment, name, tree_directory, info_directory, datatype)
