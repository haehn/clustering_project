#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import os
import dendropy as dpy
from subprocess import Popen, PIPE
from tree import Tree


class SequenceRecord(object):

    """ 
    Class for reading sequence files in fasta or phylip formats
    Supports writing in fasta, phylip, phylip interleaved, nexus formats,
    sorting sequences by length and name,
    concatenating sequences when sequence names are a perfect match,
    iterating over records.
    """

    def get_fasta_file(
        self,
        fasta_file,
        name=None,
        datatype=None,
        ):
        """ FASTA format parser: turns fasta file into Alignment_record object 
        """

        headers = []
        sequences = []
        openfile = open(fasta_file, 'r')

        # skip over file until first header is found

        while True:
            line = openfile.readline()
            if not line:
                return
            if line[0] == '>':
                break

            # we break the loop here at the start of the first record

        headers.append(line[1:].rstrip())  # chuck the first header into our list

        while True:
            line = openfile.readline()
            sequence_so_far = []  # build up sequence a line at a time in this list
            while True:
                if not line:
                    break
                elif not line[0] == '>':
                    sequence_so_far.append(line.rstrip())
                    line = openfile.readline()
                else:
                    break
            sequences.append(''.join(sequence_so_far).replace(',', ''))
            if not line:
                break
            headers.append(line[1:].rstrip())

        # check all sequences the same length

        first_seq_length = len(sequences[0])
        is_alignment = True
        for seq in sequences:
            if len(seq) != first_seq_length:
                is_alignment = False
                break
            else:
                continue

        # check same number of headers as sequences

        if len(headers) != len(sequences):
            print 'Error matching all headers and sequences'

        if is_alignment:
            self.is_aligned = True
        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.datatype = datatype
        self._update()

    def get_phylip_file(
        self,
        phylip_file,
        name=None,
        datatype=None,
        ):
        """ PHYLIP format parser"""

        mapping = {}
        openfile = open(phylip_file, 'r')
        info = openfile.readline().split()
        num_taxa = int(info[0])
        seq_length = int(info[1])

        for line in openfile:
            line = line.rstrip()
            if not line:
                continue
            line = line.split()
            header = line[0]
            sequence = line[1]
            if not header in mapping:
                mapping[header] = sequence
            else:
                mapping[header] += sequence

        assert len(mapping) == num_taxa
        for each_sequence in mapping.values():
            assert len(each_sequence) == seq_length

        self.name = name
        self.datatype = datatype
        (self.headers, self.sequences) = (mapping.keys(),
                mapping.values())
        self._update()

    def __init__(
        self,
        infile=None,
        file_format='fasta',
        name=None,
        headers=None,
        sequences=None,
        ):

        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.mapping = None
        self.length = 0
        self.seqlength = 0
        self.datatype = None
        self.is_aligned = False
        if infile:
            if file_format == 'fasta':
                self.get_fasta_file(infile, **kwargs)
            elif file_format == 'phylip':
                self.get_phylip_file(infile, **kwargs)
        self.index = -1
        self._update()

    def _update(self):
        """ For updating the length and mapping attributes of the object after 
        reading sequences """

        if self.headers and self.sequences:
            self.mapping = dict(zip(self.headers, self.sequences))
            self.length = len(self.mapping)
            first_seq_length = len(self.sequences[0])
            is_aligned = True
            for seq in self.sequences:
                if len(seq) != first_seq_length:
                    is_aligned = False
                    break
                else:
                    continue
            self.is_aligned = is_aligned
            self.seqlength = first_seq_length

    def __iter__(self):  # Should do this with generators / yield
        return self

    def next(self):  # As above
        self.index += 1
        if self.index == self.length:
            self.index = -1
            raise StopIteration
        return {'header': self.headers[self.index],
                'sequence': self.sequences[self.index]}

    def __str__(self):
        output_string = ''
        if self.is_aligned:
            output_string += \
                'Aligned_Sequence_Record: {0}\n'.format(self.name)
        else:
            output_string += \
                'Unaligned_Sequence_Record: {0}\n'.format(self.name)
        for i in range(len(self.headers)):
            if len(self.sequences[i]) > 50:
                output_string += '>{0}\n'.format(self.headers[i]) \
                    + '{0}\n'.format((self.sequences[i])[:50]) \
                    + '... ({0}) ...\n'.format(len(self.sequences[i])
                        - 100) \
                    + '{0}\n'.format((self.sequences[i])[-50:]) + '\n'
            else:
                output_string += '>{0}\n{1}'.format(self.headers[i],
                        self.sequences[i]) + '\n' + '\n'
        output_string += '{0} sequences in record'.format(len(self))
        return output_string

    def __len__(self):
        return self.length

    def __add__(self, other):
        """ 
        Allows SequenceRecords to be added together, concatenating sequences
        """

        try:
            assert set(self.headers) == set(other.headers)
        except AssertionError:
            print 'Sequence labels do not match between alignments'
            return self
        d = {}
        for k in self.mapping.keys():
            d[k] = self.mapping[k] + other.mapping[k]
        return SequenceRecord(headers=d.keys(), sequences=d.values())

    def __radd__(self, other):
        return other.__add__(self)

    def __mul__(self, n):
        """ 
        Allows SequenceRecord * x to return a concatenation of the record x 
        times
        """

        if not isinstance(n, int):
            print 'Truncating {0} to {1}'.format(n, int(n))
            n = int(n)
        new_seqs = []
        for s in self.sequences:
            new_seqs.append(s * n)
        return SequenceRecord(headers=self.headers, sequences=new_seqs)

    def __rmul__(self, n):
        if not isinstance(n, int):
            print 'Truncating {0} to {1}'.format(n, int(n))
            n = int(n)
        new_seqs = []
        for s in self.sequences:
            new_seqs.append(s * n)
        return SequenceRecord(headers=self.headers, sequences=new_seqs)

    def __eq__(self, other):
        if type(other) is type(self):
            return set(self.headers) == set(other.headers) \
                and set(self.sequences) == set(other.sequences)
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def sort_by_length(self, in_place=True):
        """
        Sorts sequences by ungapped length
        If in_place = False the sorting doesn't mutate the underlying
        object, and the output is returned
        If in_place = True the sorting mutates the self object
        """

        # Sort sequences by descending order of length
        # Uses zip as its own inverse [ zip(*zip(A,B)) == (A,B) ]

        (h, s) = zip(*sorted(zip(self.headers, self.sequences),
                     key=lambda item: len(item[1].replace('-', '')),
                     reverse=True))
        if in_place:
            self.headers = h
            self.sequences = s
        return SequenceRecord(name=self.name, headers=h, sequences=s)

    def sort_by_name(self, in_place=True):
        """
        Sorts sequences by name, treating numbers as integers (i.e. 
        sorting like this: 1, 2, 3, 10, 20 not 1, 10, 2, 20, 3).
        If in_place = False the sorting doesn't mutate the underlying
        object, and the output is returned
        If in_place = True the sorting mutates the self object
        """

        items = self.mapping.items()
        sort_key = lambda item: tuple((int(num) if num else alpha)
                for (num, alpha) in re.findall(r'(\d+)|(\D+)', item[0]))
        items = sorted(items, key=sort_key)
        (h, s) = zip(*items)
        if in_place:
            self.headers = h
            self.sequences = s
        else:
            return SequenceRecord(name=self.name, headers=h,
                                  sequences=s)

    def write_fasta(self, outfile='stdout', print_to_screen=False):
        """ 
        Writes sequences to file in fasta format
        If outfile = 'stdout' the sequences are printed to screen, not written
        to file
        If print_to_screen = True the sequences are printed to screen 
        whether they are written to disk or not
        """

        lines = ['>{0}\n{1}'.format(h, seq) for (h, seq) in
                 zip(self.headers, self.sequences)]
        s = '\n'.join(lines)
        s += '\n'
        if outfile == 'stdout':
            print s
            return s
        elif outfile == 'pipe':
            if print_to_screen:
                print s
            return s
        else:
            open(outfile, 'w').write(s)
            if print_to_screen:
                print s

    def write_nexus(self, outfile='stdout', sequence_type='protein'):
        maxlen = len(max(self.sequences, key=len))
        lines = ['{0:<14} {1:-<{2}}'.format(x, y, maxlen) for (x, y) in
                 zip(self.headers, self.sequences)]
        file_header = '#NEXUS\n' + '\n'
        file_header += 'begin data;\n'
        file_header += \
            '    dimensions ntax={0} nchar={1};\n'.format(self.length,
                maxlen)
        file_header += \
            '    format datatype={0} interleave=no gap=-;\n'.format(sequence_type)
        file_header += '    matrix\n' + '\n' + '\n'

        file_footer = '    ;\n' + '\n' + 'end;\n'

        s = file_header + '\n'.join(lines) + file_footer
        if outfile == 'stdout':
            print s
            return s
        elif outfile == 'pipe':
            return s
        else:
            open(outfile, 'w').write(s)

    def write_phylip(
        self,
        outfile='stdout',
        print_to_screen=False,
        interleaved=False,
        line_length=120,
        ):
        """ 
        Writes sequences to file in phylip format, interleaving optional
        If outfile = 'stdout' the sequences are printed to screen, not written
        to disk
        If print_to_screen = True the sequences are printed to screen 
        whether they are written to disk or not
        """

        maxlen = len(max(self.sequences, key=len))
        file_header = ' {0} {1}\n'.format(self.length, maxlen)
        if interleaved:
            maxheader = len(max(self.headers, key=len))
            num_lines = maxlen / line_length + 1
            s = file_header
            for i in range(num_lines):
                for seq_header in self.headers:
                    s += '{0:<{1}} {2}\n'.format(seq_header,
                            max(maxheader + 1, 15),
                            (self.mapping[seq_header])[i
                            * line_length:(i + 1) * line_length])
                s += '\n'
        else:
            lines = ['{0:<14} {1:-<{2}}'.format(x, y, maxlen) for (x,
                     y) in zip(self.headers, self.sequences)]
            s = file_header + '\n'.join(lines)
            s += '\n'
        if outfile == 'stdout':
            print s
            return s
        elif outfile == 'pipe':
            if print_to_screen:
                print s
            return s
        else:
            open(outfile, 'w').write(s)
            if print_to_screen:
                print s


class TCSeqRec(SequenceRecord):

    """
    A version of the SequenceRecord class with some extra functionality
    for working with tree inference packages, notably TreeCollection
    """

    def __init__(
        self,
        infile=None,
        file_format='fasta',
        name=None,
        headers=None,
        sequences=None,
        dv=[],
        datatype=None,
        tree=Tree(),
        ):

        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.datatype = datatype
        self.length = 0
        self.seqlength = 0
        self.is_aligned = False
        self.tree = tree
        self.TCfiles = {}
        self.dv = dv
        if infile:
            if file_format == 'fasta':
                self.get_fasta_file(infile, name=name,
                                    datatype=datatype)
            elif file_format == 'phylip':
                self.get_phylip_file(infile, name=name,
                        datatype=datatype)
        self.index = -1
        self._update()

    def __add__(self, other):
        """ 
        Allows Records to be added together, concatenating sequences
        """

        self_set = set(self.headers)
        other_set = set(other.headers)
        if not self.datatype == other.datatype:
            print 'Trying to add sequences of different datatypes'
            return self
        union = self_set | other_set
        intersection = self_set & other_set
        only_in_self = self_set - other_set
        only_in_other = other_set - self_set

        d = {}
        for k in union:
            if k in intersection:
                d[k] = self.mapping[k] + other.mapping[k]
            elif k in only_in_self:
                d[k] = self.mapping[k] + 'N' * other.seqlength
            elif k in only_in_other:
                d[k] = 'N' * self.seqlength + other.mapping[k]
        dvsum = self.dv + other.dv
        return_object = TCSeqRec(headers=d.keys(),
                                 sequences=d.values(),
                                 datatype=self.datatype).sort_by_name(in_place=False)
        return_object.dv = dvsum
        return return_object

    def sort_by_length(self, in_place=True):
        """
        Sorts sequences by descending order of length
        Uses zip as its own inverse [ zip(*zip(A,B)) == (A,B) ]
        Gaps and 'N' characters are not counted
        """

        (h, s) = zip(*sorted(zip(self.headers, self.sequences),
                     key=lambda item: len(item[1].replace('-', ''
                     ).replace('N', '')), reverse=True))
        if in_place:
            self.headers = h
            self.sequences = s
        else:
            return TCSeqRec(name=self.name, headers=h, sequences=s,
                            datatype=self.datatype)

    def sort_by_name(self, in_place=True):
        """
        Sorts sequences by name, treating numbers as integers (i.e. 
        sorting like this: 1, 2, 3, 10, 20 not 1, 10, 2, 20, 3).
        If in_place = False the sorting doesn't mutate the underlying
        object, and the output is returned
        If in_place = True the sorting mutates the self object
        """

        items = self.mapping.items()
        sort_key = lambda item: tuple((int(num) if num else alpha)
                for (num, alpha) in re.findall(r'(\d+)|(\D+)', item[0]))
        items = sorted(items, key=sort_key)
        (h, s) = zip(*items)
        if in_place:
            self.headers = h
            self.sequences = s
        else:
            return TCSeqRec(name=self.name, headers=h, sequences=s,
                            datatype=self.datatype)

    def sanitise(self):
        self.sort_by_name()
        l = []
        for h in self.headers:
            if '/' in h:
                h = h[:h.index('/')]
            while h.startswith(' '):
                h = h[1:]
            h = h.replace(' ', '_')
            l.append(h)
        self.headers = l
        self.sequences = [seq.upper() for seq in self.sequences]
        self._update()

    def _write_temp_phylip(self, tmpdir='/tmp'):
        self.write_phylip('{0}/{1}.phy'.format(tmpdir, self.name))

    def _write_temp_tc(self, tmpdir='/tmp'):
        num_matrices = len(self.dv)
        if num_matrices == 0:
            print 'No distance-variance matrix available'
            return
        all_labels = self.headers
        len_all_labels = len(all_labels)

        # Temporary files

        dv_tmpfile = open('{0}/{1}_dv.txt'.format(tmpdir, self.name),
                          'w')
        labels_tmpfile = open('{0}/{1}_labels.txt'.format(tmpdir,
                              self.name), 'w')
        map_tmpfile = open('{0}/{1}_map.txt'.format(tmpdir, self.name),
                           'w')
        guidetree_tmpfile = open('{0}/{1}_tree.nwk'.format(tmpdir,
                                 self.name), 'w')

        # Write headers to temp files

        dv_tmpfile.write('{0}\n'.format(num_matrices))
        map_tmpfile.write('{0} {1}\n'.format(num_matrices,
                          len_all_labels))
        labels_tmpfile.write('''{0}
{1}
'''.format(len_all_labels,
                             ' '.join(all_labels)))
        labels_tmpfile.flush()

        # Add info from self.dv to temp files

        for (i, (matrix, labels)) in enumerate(self.dv):
            labels = labels.split()
            dim = len(labels)
            index = i + 1
            dv_tmpfile.write('''{0} {0} {1}
{2}
'''.format(dim, index,
                             matrix))
            dv_tmpfile.flush()
            for lab in all_labels:
                if lab in labels:
                    map_tmpfile.write('{0} '.format(labels.index(lab)
                            + 1))
                else:
                    map_tmpfile.write('-1 ')
            map_tmpfile.write('\n')
            map_tmpfile.flush()

        # Close finished temp files

        map_tmpfile.close()
        dv_tmpfile.close()
        labels_tmpfile.close()

        # Write the guidetree

        guidetree = self.get_guide_tree(tmpdir)
        guidetree = guidetree.newick
        guidetree_tmpfile.write(guidetree)
        guidetree_tmpfile.flush()
        guidetree_tmpfile.close()

        # Check it all worked

        assert os.path.isfile('{0}/{1}_dv.txt'.format(tmpdir, self.name))
        assert os.path.isfile('{0}/{1}_labels.txt'.format(tmpdir, self.name))
        assert os.path.isfile('{0}/{1}_map.txt'.format(tmpdir, self.name))
        assert os.path.isfile('{0}/{1}_tree.nwk'.format(tmpdir, self.name))

    def get_phyml_tree(self, tmpdir='/tmp'):
        self.tree = Tree()
        self._write_temp_phylip(tmpdir=tmpdir)
        print 'Running phyml on ' + str(self.name) + '...'
        input_file = '{0}/{1}.phy'.format(tmpdir, self.name)
        if self.datatype == 'dna':
            model = 'GTR'
            datatype = 'nt'
        elif self.datatype == 'protein':
            model = 'WAG'
            datatype = 'aa'
        else:
            print 'I don\'t know this datatype: {0}'.format(self.datatype)
            return
        t = self.tree.run_phyml(model, input_file, datatype, self.name)
        os.remove('{0}/{1}.phy'.format(tmpdir, self.name))
        return self.tree

    def get_raxml_tree(self, tmpdir='/tmp'):
        self.tree = Tree()
        self._write_temp_phylip(tmpdir=tmpdir)
        print 'Running raxml on ' + str(self.name) + '...'
        input_file = '{0}/{1}.phy'.format(tmpdir, self.name)
        if self.datatype == 'dna':
            model = 'GTRGAMMA'
        elif self.datatype == 'protein':
            model = 'PROTGAMMAWAG'
        else:
            print 'I don\'t know this datatype: {0}'.format(self.datatype)
            return
        self.tree.run_raxml(model, input_file, self.name, tmpdir)
        os.remove('{0}/{1}.phy'.format(tmpdir, self.name))
        if os.path.isfile('{0}/{1}.phy.reduced'.format(tmpdir,
                          self.name)):
            os.remove('{0}/{1}.phy.reduced'.format(tmpdir, self.name))
        return self.tree

    def get_guide_tree(self, tmpdir='/tmp'):
        self._write_temp_phylip(tmpdir=tmpdir)
        input_file = '{0}/{1}.phy'.format(tmpdir, self.name)
        if self.datatype == 'dna':
            model = 'GTRGAMMA'
        elif self.datatype == 'protein':
            model = 'PROTGAMMAWAG'
        else:
            print 'I don\'t know this datatype: {0}'.format(self.datatype)
            return
        t = Tree().run_raxml(model, input_file, self.name, tmpdir,
                             guide=True)
        if os.path.isfile('{0}/{1}.phy'.format(tmpdir, self.name)):
            os.remove('{0}/{1}.phy'.format(tmpdir, self.name))
        if os.path.isfile('{0}/{1}.phy.reduced'.format(tmpdir,
                          self.name)):
            os.remove('{0}/{1}.phy.reduced'.format(tmpdir, self.name))
        return t

    def get_TC_tree(self, tmpdir='/tmp'):
        self._write_temp_tc(tmpdir=tmpdir)
        print 'Running TreeCollection on ' + str(self.name) + '...'
        self.tree = \
            Tree().run_treecollection('{0}/{1}_dv.txt'.format(tmpdir,
                self.name), '{0}/{1}_map.txt'.format(tmpdir,
                self.name), '{0}/{1}_labels.txt'.format(tmpdir,
                self.name), '{0}/{1}_tree.nwk'.format(tmpdir,
                self.name), self.name)
        os.remove('{0}/{1}_dv.txt'.format(tmpdir, self.name))
        os.remove('{0}/{1}_map.txt'.format(tmpdir, self.name))
        os.remove('{0}/{1}_labels.txt'.format(tmpdir, self.name))
        os.remove('{0}/{1}_tree.nwk'.format(tmpdir, self.name))
        return self.tree

    def get_dv_matrix(self, tmpdir='/tmp',
                      helper='/Users/kgori/Projects/clustering_project/class_files/DV_wrapper.drw'
                      ):
        """
        Makes a call to the TC_wrapper.drw darwin helper script, which
        calculates a distance-variance matrix from the sequence alignments,
        and generates files needed by the treecollection binary
        """

        if self.name:
            fastafile = '{0}/{1}.fas'.format(tmpdir, self.name)
        else:
            fastafile = '{0}/fasta_tmp.fas'.format(tmpdir)
        if self.datatype == 'dna':
            datatype = 'DNA'
        else:
            datatype = 'AA'

        self.write_fasta(fastafile)
        command = \
            'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\'; fpath := \'{2}/\'; ReadProgram(\'{3}\');" | darwin'.format(fastafile,
                datatype, tmpdir, helper)
        process = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
        (stdout, stderr) = process.communicate()
        dv_string = \
            open('{0}/temp_distvar.txt'.format(tmpdir)).read().rstrip()
        labels = ' '.join(self.headers)
        os.remove(fastafile)
        os.remove('{0}/temp_distvar.txt'.format(tmpdir))
        return (dv_string, labels)
