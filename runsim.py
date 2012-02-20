#!/usr/bin/env python
"""
Quick & dirty script to run alfsim and generate 25 genes for 20 species, with genes
from two classes (i.e. a mixture of two defined histories - 15 from class 1 and 10 
from class 2).
MSAs are put in MSA directory
True trees are put in trees directory
Remaining ALF files are deleted
"""
from ALF_wrapper import *         # functions write_ALF_parameters() and run_ALF()
from sequence_record import *     # classes Unaligned_Sequence_Record and Aligned_Sequence_Record, function get_fasta_file()
from handleArgs import handleArgs
import glob, os, shutil, sys

#Write parameters
write_ALF_parameters('class1', 'scratch', '', 15, 250, 20, 250, True, 'class1-params.drw')
write_ALF_parameters('class2', 'scratch', '', 10, 250, 20, 250, True, 'class2-params.drw')

#Run simulation
run_ALF('class1-params.drw')
run_ALF('class2-params.drw')

#Gather MSA files
MSA_path = "./MSA"

if not os.path.isdir(MSA_path): os.mkdir(MSA_path)
class1_path = './scratch/class1/MSA'
class2_path = './scratch/class2/MSA'

class1_msas = glob.glob("{0}/*.fa".format(class1_path))
class1_msas = sorted(class1_msas, key=lambda name: int(name.split('_')[1]))

class2_msas = glob.glob("{0}/*.fa".format(class2_path))
class2_msas = sorted(class2_msas, key=lambda name: int(name.split('_')[1]))

for fi in class1_msas:
    name_elements = fi.split('_')
    rename = "{0}/MSA_{1:0>2}.fa".format(MSA_path, int(name_elements[1]))
    print fi, rename
    os.rename(fi, rename)

for fi in class2_msas:
    name_elements = fi.split('_')
    rename = "{0}/MSA_{1:0>2}.fa".format(MSA_path, int(name_elements[1])+15)
    print fi, rename
    os.rename(fi, rename)

#Gather tree files
tree_path = "./trees"
if not os.path.isdir(tree_path): os.mkdir(tree_path)
os.rename('./scratch/class1/RealTree.nwk', './trees/true1.nwk')
os.rename('./scratch/class2/RealTree.nwk', './trees/true2.nwk')

#Delete unnecessary simulated filesls
shutil.rmtree("./scratch/class1")
shutil.rmtree("./scratch/class2")

#Convert fasta files to phylip files for use with raxml
for fasta in glob.glob("{0}/*.fa".format(MSA_path)):
    phylip = fasta[:fasta.rindex('.')]+'.phy'
    seq_object = get_fasta_file(fasta)
    seq_object.headers = [x.split('/')[0] for x in seq_object.headers]
    seq_object.write_phylip(outfile=phylip, print_to_screen = True)
