#!/usr/bin/env python

from Bio import Entrez
from pprint import pprint
Entrez.email = 'kgori@ebi.ac.uk'

inf = open('gene_table.tsv')
data = [x.rstrip().split("\t") for x in inf.readlines()]

labels = data[0]
data = data[1:]

nrows = len(data)
ncols = len(labels)

loci = {}
sequences = {}

for i in range(2,ncols):
    for j in range(2,nrows):
        if data[j][i] not in ['','XXXXXX','XXXXXXX','n/a']:
            if labels[i] not in loci: loci[labels[i]] = [data[j][i]]
            else: loci[labels[i]].append(data[j][i])
            if data[j][i] not in sequences:
                taxon = '_'.join(data[j][0].split()[::len(data[j][0].split())-1])
                label = data[j][1].replace(' ','_')
                print "Fetching record {0} ...".format(data[j][i])
                sequence = '\n'.join(Entrez.efetch(db="nucleotide", id=data[j][i], rettype="fasta", retmode="text").read().split('\n')[1:])
                sequences[data[j][i]] = {'Taxon': taxon, 'Label': label, 'Sequence': sequence}

for locus in loci:
    writer = open ("Unaligned_Sequences/{0}.fas".format(locus),"w")
    for seq_id in loci[locus]:
        print "writing record {0},{1} ...".format(locus,seq_id)
        taxon = sequences[seq_id]['Taxon']
        label = sequences[seq_id]['Label']
        sequence = sequences[seq_id]['Sequence']
        if label != 'n/a': output = ">{0}_({1})\n{2}".format(taxon,label,sequence)
        else: output = ">{0}\n{1}".format(taxon,sequence)
        writer.write (output)
    writer.close()