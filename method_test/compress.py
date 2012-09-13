#!/usr/bin/env python
import argparse
import os
import sys
import glob
import tarfile

def fpath(s):
    """
    Helper function used when passing filepath arguments with argparse module.
        Trims all '/' characters from the end of the path string.
    """
    while s.endswith('/'):
        s = s[:-1]
    return s

def getFolderSize(folder):
    total_size = os.path.getsize(folder)
    for item in os.listdir(folder):
        itempath = os.path.join(folder, item)
        if os.path.isfile(itempath):
            total_size += os.path.getsize(itempath)
        elif os.path.isdir(itempath):
            total_size += getFolderSize(itempath)
    return total_size

def sizemessage(indir):
   print '{0} takes up {1} bytes'.format(indir,getFolderSize(indir))

sort_key = lambda item: tuple((int(num) if num else alpha) for (num,alpha) in re.findall(r'(\d+)|(\D+)', item))

parser = argparse.ArgumentParser(prog='compress.py')
parser.add_argument('-d', '--directory', help='input directory', type=fpath, default='.')
args = vars(parser.parse_args())
indir = args['directory']

dna_dir = '{0}/dna_alignments'.format(indir)
aa_dir = '{0}/aa_alignments'.format(indir)
bionj_dir = '{0}/bionj_clustering'.format(indir)
phyml_dir = '{0}/phyml_clustering'.format(indir)

dna_fasta_files = glob.glob('{0}/*.fas'.format(dna_dir))
aa_fasta_files = glob.glob('{0}/*.fas'.format(aa_dir))
dna_phylip_files = glob.glob('{0}/*.phy'.format(dna_dir))
aa_phylip_files = glob.glob('{0}/*.phy'.format(aa_dir))
bionj_phylip_files = glob.glob('{0}/*.phy'.format(bionj_dir))
phyml_phylip_files = glob.glob('{0}/*.phy'.format(phyml_dir))



sizemessage(indir)
print 'Deleting fasta files...'
for fasta_file in dna_fasta_files + aa_fasta_files:
    os.remove(fasta_file)

sizemessage(indir)
print 'Compressing dna alignments...'
with tarfile.open('{0}/dna_alignments.tar.bz2'.format(dna_dir), 'w:bz2') as dnatar:
    for phylip_file in dna_phylip_files:
        dnatar.add(phylip_file)
        os.remove(phylip_file)

sizemessage(indir)
print 'Compressing aa alignments...'
with tarfile.open('{0}/aa_alignments.tar.bz2'.format(aa_dir), 'w:bz2') as aatar:
    for phylip_file in aa_phylip_files:
        aatar.add(phylip_file)
        os.remove(phylip_file)

sizemessage(indir)
print 'Compressing bionj clustering alignments...'
with tarfile.open('{0}/bionj_clustering_alignments.tar.bz2'.format(bionj_dir), 'w:bz2') as bionjtar:
    for phylip_file in bionj_phylip_files:
        bionjtar.add(phylip_file)
        os.remove(phylip_file)

sizemessage(indir)
print 'Compressing phyml clustering alignments...'
with tarfile.open('{0}/phyml_clustering_alignments.tar.bz2'.format(phyml_dir), 'w:bz2') as phymltar:
    for phylip_file in phyml_phylip_files:
        phymltar.add(phylip_file)
        os.remove(phylip_file)
sizemessage(indir)

