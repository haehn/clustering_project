#!/usr/bin/env python
"""
Script to run alfsim and generate G genes for 20 species, with genes from C classes (i.e. a mixture of defined histories).
MSAs are put in MSA directory
True trees are put in trees directory
Remaining ALF files are deleted
"""
# from ALF_wrapper import *         # functions write_ALF_parameters() and run_ALF()
from sequence_record import *     # classes Unaligned_Sequence_Record and Aligned_Sequence_Record, function get_fasta_file()
# from utility_functions import * 
import argparse, glob, os, re, shutil, sys

####################################################################################################

class Unaligned_Sequence_Record(object):
    """ Class for holding unaligned sequences """
    def __init__(self, name, headers, sequences):
        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.mapping = dict(zip(headers,sequences)) # Probably can lose this
        self.index = -1
        self.length = len(self.headers)
    def __iter__(self): # Should do this with generators / yield
        return self
    def next(self): # As above
        self.index += 1
        if self.index == self.length: 
            self.index = -1
            raise StopIteration
        return { 'header' : self.headers[self.index],
                 'sequence' : self.sequences[self.index]}
    def __str__(self):
        output_string = ''
        output_string += 'Unaligned_Sequence_Record: {0}\n'.format( self.name )
        for i in range(len(self.headers)):
            output_string += '>{0}\n{1}...({2})'.format( self.headers[i], self.sequences[i][:50],len(self.sequences[i]) ) + '\n'
        output_string += "{0} sequences in record".format(len(self))
        return output_string
    def __len__(self):
        return self.length

    def sort_by_length(self):
        # Sort sequences by descending order of length
        # Uses zip as its own inverse [ zip(*zip(A,B)) == (A,B) ]
        h, s = zip(*sorted(zip(self.headers,self.sequences), key = lambda item: len(item[1]),reverse=True))
        self.headers = h
        self.sequences = s
        return self

    def sort_by_name(self):
        h, s = zip(*sorted(zip(self.headers,self.sequences)))
        self.headers = h
        self.sequences = s
        return self

    def write_fasta(self,outfile="stdout",print_to_screen=False):
        lines = ['>{0}\n{1}'.format(h,seq) for h,seq in zip(self.headers,self.sequences)]
        s = '\n'.join(lines)
        s += '\n'
        if outfile == "stdout":
            print s
            return s
        elif outfile == "pipe":
            if print_to_screen: print s
            return s
        else: 
            open(outfile,'w').write(s)
            if print_to_screen: print s

    def write_nexus(self,outfile="stdout",sequence_type='protein'):
        maxlen = len(max(self.sequences,key=len))
        lines = [ "{0:<14} {1:-<{2}}".format(x,y,maxlen) for (x,y) in zip(self.headers,self.sequences) ]
        file_header =  "#NEXUS\n\n"
        file_header += "begin data;\n"
        file_header += "    dimensions ntax={0} nchar={1};\n".format(self.length,maxlen)
        file_header += "    format datatype={0} interleave=no gap=-;\n".format(sequence_type)
        file_header += "    matrix\n\n"

        file_footer = "    ;\n\nend;\n"

        s = file_header + '\n'.join(lines) + file_footer
        if outfile == "stdout":
            print s
            return s
        elif outfile == "pipe":
            return s
        else: open(outfile,'w').write(s)

    def write_phylip(self,outfile="stdout",sequence_type='protein',print_to_screen=False):
        maxlen = len(max(self.sequences,key=len))
        lines = [ "{0:<14} {1:-<{2}}".format(x,y,maxlen) for (x,y) in zip(self.headers,self.sequences) ]
        file_header = " {0} {1}\n".format(self.length,maxlen)
        s = file_header + '\n'.join(lines)
        s += '\n'
        if outfile == "stdout":
            print s
            return s
        elif outfile == "pipe":
            if print_to_screen: print s
            return s
        else: 
            open(outfile,'w').write(s)
            if print_to_screen: print s


class Aligned_Sequence_Record(Unaligned_Sequence_Record):
    """ Class for holding a sequence alignment """
    def __str__(self):
        output_string = ""
        output_string += "Aligned_Sequence_Record: {0}\n".format( self.name )
        for i in range(len(self.headers)):
            output_string += ">{0}\n{1}...({2})\n".format( self.headers[i], self.sequences[i][:50], len(self.sequences[i]) )
        output_string += "{0} sequences in record".format(len(self))
        return output_string
 
    def sort_by_length(self):
        h, s = zip(*sorted(zip(self.headers,self.sequences), key = lambda item: len(item[1].replace('-','')),reverse=True))
        self.headers = h
        self.sequences = s
        return self

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

def pam2sps(tree_file, conversion, outfile=None):
    """ Take in a newick file, change branch lengths by factor 100, write output file
        (i.e. convert between PAM units and substitutions per site)
    """
    import re
    reg_ex = re.compile('(?<=:)[0-9.]+')
    convert_pam_to_sps = lambda a: str(0.01*float(a.group()))
    convert_sps_to_pam = lambda b: str(100*float(b.group()))

    input_string = open(tree_file).read()
    if conversion == 'pam2sps':
        output_string = reg_ex.sub(convert_pam_to_sps,input_string)
    elif conversion == 'sps2pam':
        output_string = reg_ex.sub(convert_sps_to_pam,input_string)
    else: output_string = reg_ex.sub(convert_pam_to_sps,input_string)

    if outfile:
        writer = open(outfile,'w')
        writer.write(output_string)
        writer.close()
    else: print output_string

def write_ALF_parameters(simulation_name, experiment_directory, 
    simulation_directory,  number_of_genes=20, min_gene_length=100, 
    number_of_species=20, mutation_rate=200, indels=True, bd_tree=True,
    gtr=False, custom_tree=False, output_filename="alf-params.drw"):
    """ Function to write parameter files for ALF """
    import os
    alfsim_parameter_string = '''# This parameter file has been generated by the ALF web service.
# To run the simulation on your own machine, download ALF from http://www.cbrg.ethz.ch/alf
# and call alfsim from the parent directory of this file.

webRequest := false;
uuid := '';

# name of simulation - you may want to change this
mname := {0};

# directories for file storage - you may want to change these
wdir := '{1}/{2}';
dbdir := 'DB/';
dbAncdir := 'DBancestral/';

# time scale for simulation (PAM is default)
unitIsPam := true:

# parameters concerning the root genome
realseed := false;
protStart := {3};
minGeneLength := {4};
gammaLengthDist := [3, 133.8063];
blocksize := 1:

# parameters concerning the substitution models
substModels := [SubstitutionModel('WAG')];
indelModels := [IndelModel(0.00001, ZIPF, [1.821], 10)];
rateVarModels := [RateVarModel()];
modelAssignments := [1]:
modelSwitchS := [[1]]:
modelSwitchD := [[1]]: '''.format( simulation_name, experiment_directory, simulation_directory, number_of_genes, min_gene_length)
    bd_tree_string = '''

# parameters concerning the species tree
treeType := 'BDTree';
birthRate := 0.01;
deathRate := 0.001;
NSpecies := {0};
ultrametric := false;
mutRate := {1};
scaleTree := false; '''.format( number_of_species, mutation_rate )
    user_tree_string = '''
# parameters concerning the species tree
treeType := 'Custom';
treeFile := '{0}'; '''.format( custom_tree )

    if not indels:
        alfsim_parameter_string = alfsim_parameter_string.replace("indelModels := [IndelModel(0.0001, ZIPF, [1.821], 10)];\n","")
    if gtr: 
        alfsim_parameter_string = alfsim_parameter_string.replace("[SubstitutionModel('WAG')];","substModels := [SubstitutionModel('GTR', [14.9, 1.42, 1, 3.04, 2.44, 8.46], [0.31,0.18,0.21,0.30], false)];")
    if bd_tree:
        alfsim_parameter_string += bd_tree_string
    if custom_tree:
        alfsim_parameter_string += user_tree_string
    if output_filename:
        if os.path.isfile(output_filename):
            write=raw_input("Output file '{0}' exists, overwrite (y/n)?: ".format(output_filename))
        else: write = 'y'
        if write == 'y':
            output = open(output_filename, 'w')
            output.write(alfsim_parameter_string)
            output.close()
    return alfsim_parameter_string


def run_ALF(parameters,quiet=False):
    """ Function to run ALF from parameter file """
    import os
    
    if os.path.isfile(parameters):
        if quiet: 
            os.system("alfsim {0} > /dev/null 2> /dev/null".format(parameters))
        else: 
            os.system("alfsim {0}".format(parameters))
    else: print "Can't find file '{0}'".format(parameters)

def get_fasta_file(fasta_file, name = "no name"):
    """ FASTA format parser: turns fasta file into Alignment_record object """
    headers = []
    sequences = []
    openfile = open(fasta_file,'r')

    #skip over file until first header is found
    while True:
        line = openfile.readline()
        if not line: return
        if line[0] == ">": break
        #we break the loop here at the start of the first record
    
    headers.append(line[1:].rstrip()) #chuck the first header into our list

    while True:
        line = openfile.readline()
        sequence_so_far = [] #build up sequence a line at a time in this list
        while True:
            if not line: break
            elif not line[0] == ">": 
                sequence_so_far.append(line.rstrip())
                line = openfile.readline()    
            else: break
        sequences.append("".join(sequence_so_far).replace(",",""))
        if not line: break
        headers.append(line[1:].rstrip())
    
    #check all sequences the same length
    first_seq_length = len(sequences[0])
    is_alignment = True
    for seq in sequences:
        if len(seq) != first_seq_length:
            is_alignment = False
            break
        else: continue

    #check same number of headers as sequences
    if len(headers) != len(sequences): print "Error matching all headers and sequences"

    if is_alignment: return Aligned_Sequence_Record(name, headers, sequences)
    else: return Unaligned_Sequence_Record(name, headers, sequences)


#=============================#
# get command line arguments #
#===========================#

parser = argparse.ArgumentParser(prog='runsim.py', description='Run ALF simulator to generate different topological classes')
parser.add_argument('-c','--classes', help='Number of classes', type=int, default=2)
parser.add_argument('-s','--species', help='Number of species', type=int, default=20)
parser.add_argument('-g','--genes', help='Number of genes per class', nargs='*', default=[])
parser.add_argument('-R','--REV', help="Use GTR / REV model", dest='gtr', action='store_true')
parser.add_argument('-t','--custom_trees', dest='custom_trees', help="Use custom tree(s)", nargs='*', default=[])
parser.add_argument('-P','--PAM',help='Convert sps branch lengths to PAM',action='store_true',dest='pam')
parser.add_argument('-S','--SPS',help='Convert PAM branch lengths to sps',action='store_true',dest='sps')
parser.add_argument('-r','--rates', help='Tree length per class (PAM units)', nargs='*', default=[])
parser.add_argument('-d','-dir','--directory', help='Base output directory', type=fpath, default='.')
parser.add_argument('-tmp','--temp-directory', help='Directory to use for temp files', type=fpath, default='./alftmp')
parser.add_argument('-q','--quiet',dest='quiet',action='store_true',help='Not much printing to stdout')
args = vars(parser.parse_args())

# Do some tidying up
while len(args['genes']) < args['classes']: # Make sure each class has a number of genes associated with it
    args['genes'].append(10)
while len(args['rates']) < args['classes']: # Make sure each class has a rate
    args['rates'].append(200)

# Shorter variable names (CAPS = command line args)
C = args['classes']
S = args['species']
G = [int(x) for x in args['genes']]
R = [int(y) for y in args['rates']]
OUT_DIR = args['directory']
TEMP_DIR = args['temp_directory']
MSA_path = fpath("{0}/MSA".format(OUT_DIR))
tree_path = fpath("{0}/trees".format(OUT_DIR))
custom_trees = args['custom_trees']
quiet = args['quiet']
gtr = args['gtr']
pam = args['pam']
sps = args['sps']
for each in [OUT_DIR, TEMP_DIR, MSA_path, tree_path]:
    if not os.path.isdir(each): os.mkdir(each)

####################################################################################################

#============#
# MAIN LOOP #
#==========#

for i in range(C):
    #Write parameters
    if custom_trees:
        try: 
            custom_tree = custom_trees[i]
        except IndexError:
            custom_tree = None

    if custom_tree:  
        if not os.path.isfile(custom_tree):
            print 'Couldn\'t open file {0}'.format(custom_tree)
            sys.exit()  
        if sps:
            pam2sps(custom_tree,'pam2sps','{0}_sps.nwk'.format(custom_tree[:custom_tree.rindex('.')]))
            custom_tree = '{0}_sps.nwk'.format(custom_tree[:custom_tree.rindex('.')])
        elif pam:
            pam2sps(custom_tree,'sps2pam','{0}_pam.nwk'.format(custom_tree[:custom_tree.rindex('.')]))
            custom_tree = '{0}_pam.nwk'.format(custom_tree[:custom_tree.rindex('.')])
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR,\
            simulation_directory='', number_of_genes=G[i], min_gene_length=130, number_of_species=S,\
            mutation_rate=R[i], indels=True, gtr=gtr, bd_tree=False, custom_tree=custom_tree,\
            output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))
    else:
        write_ALF_parameters(simulation_name='class{0}'.format(i+1), experiment_directory=TEMP_DIR,\
            simulation_directory='', number_of_genes=G[i], min_gene_length=130, number_of_species=S,\
            mutation_rate=R[i], indels=True, gtr=gtr, output_filename='{0}/class{1}-params.drw'.format(OUT_DIR,i+1))   

    #Run simulation
    run_ALF('{0}/class{1}-params.drw'.format(OUT_DIR,i+1),quiet)

    #Gather MSA files
    class_path = './{0}/class{1}/MSA'.format(TEMP_DIR,i+1)
    if not gtr:
        class_msas = glob.glob("{0}/*.fa".format(class_path))
    else:
        class_msas = glob.glob("{0}/*dna.fa".format(class_path))
    class_msas = sorted(class_msas, key=lambda name: int(name.split('_')[1]))
    for fi in class_msas: # Give MSA file a better name
        name_elements = fi.split('_')
        gene_number = int(name_elements[1]) + sum(G[:i]) #add on number of genes already named in previous classes
        rename = "{0}/gene{1:0>3}.fas".format(MSA_path, gene_number)
        print fi, rename
        seq_object = get_fasta_file(fi)
        seq_object.headers = [x.split('/')[0] for x in seq_object.headers]
        seq_object.write_fasta(rename)
        os.remove(fi)

    #Gather tree files
    pam2sps('./{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1),"pam2sps",'./{0}/true{1}_sps.nwk'.format(tree_path,i+1))
    os.rename('./{0}/class{1}/RealTree.nwk'.format(TEMP_DIR,i+1), './{0}/true{1}.nwk'.format(tree_path,i+1))
        
    #Delete temporary files
    shutil.rmtree("./{0}/".format(TEMP_DIR))

####################################################################################################

true_clustering = open("./{0}/true_clustering".format(MSA_path),"w")
for i in range(len(G)):
    for j in range(G[i]):
        true_clustering.write(str(i+1)+" ")
true_clustering.close()

