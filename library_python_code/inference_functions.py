#!/usr/bin/env python

class Inference_Result(object):
    """ Class for storing the results of phylogenetic inference """ 
    def __init__(self, tree=None, score=None, program=None, name=None):
        self.tree = tree
        self.score = score
        self.program = program
        self.name = name

    def __str__(self):
        s = "[Inference Result Object:\n"
        if self.name: s += "Name:\t" + self.name + "\n"
        s += "Program:\t{0}\nScore:\t{1}\nTree:\t]{2}\n".format(self.program,self.score,self.tree)
        return s
   
    def pam2sps(self, conversion="sps2pam"):
        import re
        reg_ex = re.compile('(?<=:)[0-9.]+')
        
        convert_pam_to_sps = lambda a: str(0.01*float(a.group()))
        convert_sps_to_pam = lambda b: str(100*float(b.group()))
        strip_lengths      = lambda c: ''

        input_string = self.tree
        if conversion == 'pam2sps':
            output_string = reg_ex.sub(convert_pam_to_sps,input_string)
        elif conversion == 'sps2pam':
            output_string = reg_ex.sub(convert_sps_to_pam,input_string)
        else: output_string = reg_ex.sub(strip_lengths,input_string).replace(':','')

        return Inference_Result(output_string, self.score, self.program, self.name)

    def read_from_file(self, infile, name=None):
        program = None
        tree = None
        score = None
        self.name = name
        reader = open(infile)
        try:
            for line in reader:
                line = [l.rstrip().replace("]","") for l in line.split()]
                if not name and line[0] == "Name:": self.name = line[1]
                elif line[0] == "Program:": self.program = line[1]
                elif line[0] == "Tree:": self.tree = line[1]
                elif line[0] == "Score:": self.score = line[1]
        except IndexError: return
            
    def write_to_file(self, outfile):
        writer = open(outfile,'w')
        writer.write(str(self))
        writer.close()

###################################################################################################

def run_raxml(model, alignment_file, name, binary="raxml", nthreads=2, verbose=False): 
    import re, os
    command = '{0} -T {1} -m {2} -s {3} -n {4} -p 121 > /dev/null '.format(binary, nthreads, model, alignment_file, name)
    if verbose: command = command.replace("> /dev/null","")
    os.system( command )
    tree = open("RAxML_bestTree.{0}".format(name)).read()
    score = float(re.compile( "(?<=Score of best tree ).+" ).search( open("RAxML_info.{0}".format(name)).read() ).group())
    os.system( "rm *.{0}".format(name) ) # Cleanup
    return Inference_Result(tree,score,binary,name)

def run_phyml(model, alignment_file, name, datatype, binary="phyml", verbose=False):
    import re, os
    command = '{0} -m {1} -i {2} -d {3} -a e -b 0 --sequential > /dev/null'.format(binary, model, alignment_file, datatype)
    if verbose: command = command.replace("> /dev/null","")
    os.system( command )
    tree = open("{0}_phyml_tree.txt".format(alignment_file)).read()
    score = float(re.compile( "(?<=Log-likelihood: ).+" ).search( open("{0}_phyml_stats.txt".format(alignment_file)).read() ).group())
    os.system( "rm {0}_phyml_tree.txt {0}_phyml_stats.txt".format(alignment_file) ) # Cleanup
    return Inference_Result(tree,score,binary,name)

def run_treecollection_from_fasta(alignment_file, name, datatype, helper="./TC_wrapper.drw", binary="TreeCollection", verbose=False):
    import os, subprocess
    # Requires darwin helper script to generate temporary TC input files from fasta alignment
    try: assert os.path.isfile(helper)
    except AssertionError: 
        print "Can't open TC_wrapper.drw"
        return 
    command = 'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\' ; ReadProgram(\'{2}\');" | darwin > /dev/null'.format(alignment_file, datatype, helper)
    if verbose: command = command.replace("> /dev/null","")
    os.system( command )

    # Input files temp_* should exist thanks to helper script
    TC_command = '{0} temp_distvar.txt temp_map.txt temp_labels.txt temp_tree.nwk'.format(binary)
    process = subprocess.Popen( TC_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    stdout, stderr = process.communicate()
    info = stdout.split()
    tree = info[-2]
    score = info[-1]
    os.system( 'rm temp_*') # Cleanup
    return_object = Inference_Result(tree,score,binary,name).pam2sps("pam2sps")
    return return_object

def run_treecollection(distvar_file, map_file, label_file, guide_tree_file, name, binary="TreeCollection", verbose=False):
    # Requires that TC input files already exist
    import subprocess
    command = '{0} {1} {2} {3} {4}'.format(binary, distvar_file, map_file, label_file, guide_tree_file)
    process = subprocess.Popen( command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    stdout, stderr = process.communicate()
    info = stdout.split()
    tree = info[-2]
    score = info[-1]
    return_object = Inference_Result(tree,score,binary,name).pam2sps("pam2sps")
    if verbose: print command, stdout
    return return_object

##################################################################################################
