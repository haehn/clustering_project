#!/usr/bin/env python
import_debugging = False
if import_debugging: print 'sequence_record imports:'
import re
if import_debugging: print '  re (sr)'
import os
if import_debugging: print '  os (sr)'
import dendropy as dpy
if import_debugging: print '  dendropy (sr)'
from subprocess import Popen, PIPE
if import_debugging: print '  subprocess::Popen, PIPE (sr)'
from tree import Tree
if import_debugging: print '  tree::Tree (sr)'
import hashlib
if import_debugging: print '  hashlib (sr)'
from random import shuffle as shf

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

        headers.append(line[1:].rstrip())  # chuck the first header into list

        while True:
            line = openfile.readline()
            sequence_so_far = []  # build up sequence a line at a time in list
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

        openfile = open(phylip_file, 'r')
        info = openfile.readline().split()
        num_taxa = int(info[0])
        seq_length = int(info[1])

        # Initialise lists to hold the headers and sequences

        headers = [None] * num_taxa
        sequences = [None] * num_taxa

        i = 0  # counter monitors how many (informative) lines we've seen
        for line in openfile:
            line = line.rstrip()
            if not line:
                continue  # skip empty lines and don't increment counter

            line = line.split()

            # IF header is None, split line into header / sequence pair
            # ELSE header is not None, we've already seen it
            # so this file must be interleaved, so carry on
            # adding sequence fragments

            if not headers[i % num_taxa]:
                header = line[0]
                sequence_fragment = ''.join(line[1:])

                headers[i % num_taxa] = header
                sequences[i % num_taxa] = [sequence_fragment]
            else:

                sequence_fragment = ''.join(line)
                sequences[i % num_taxa].append(sequence_fragment)

            i += 1  # increment counter

        sequences = [''.join(x) for x in sequences]

        # checks

        try:
            assert len(headers) == len(sequences) == num_taxa
            for sequence in sequences:
                assert len(sequence) == seq_length
        except AssertionError:
            print 'Error reading file'
            return

        self.name = name
        self.datatype = datatype
        (self.headers, self.sequences) = (headers, sequences)
        self._update()

    def __init__(
        self,
        infile=None,
        file_format='fasta',
        name=None,
        datatype=None,
        headers=[],
        sequences=[],
        ):

        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.mapping = {}
        self.length = 0
        self.seqlength = 0
        self.datatype = datatype
        self.is_aligned = False
        if infile:
            if file_format == 'fasta':
                self.get_fasta_file(infile, name=name,
                                    datatype=datatype)
            elif file_format == 'phylip':
                self.get_phylip_file(infile, name=name,
                        datatype=datatype)
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

    @staticmethod
    def linebreaker(string, length):
        for i in range(0, len(string), length):
            yield string[i:i + length]

    def make_chunks(self, chunksize):
        num_chunks = self.seqlength / chunksize
        if self.seqlength % chunksize:
            num_chunks += 1

        new_records = []
        generators = [self.linebreaker(s, chunksize) for s in
                      self.sequences]
        for _ in range(num_chunks):
            new_record = type(self)(headers=self.headers)
            new_record.sequences = [next(g) for g in generators]
            new_record._update()
            new_records.append(new_record)
        return new_records

    def write_fasta(
        self,
        outfile='stdout',
        print_to_screen=False,
        linebreaks=None,
        ):
        """
        Writes sequences to file in fasta format
        If outfile = 'stdout' the sequences are printed to screen, not written
        to file
        If print_to_screen = True the sequences are printed to screen
        whether they are written to disk or not
        """

        if linebreaks:
            try:
                int(linebreaks)
            except ValueError:
                print 'Can\'t use {0} as value for linebreaks'.format(linebreaks)
            sequences = ['\n'.join(self.linebreaker(s, linebreaks))
                         for s in self.sequences]
        else:
            sequences = self.sequences

        lines = ['>{0}\n{1}'.format(h, seq) for (h, seq) in
                 zip(self.headers, sequences)]
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
            return outfile

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
            return outfile

    def write_phylip(
        self,
        outfile='stdout',
        print_to_screen=False,
        interleaved=False,
        linebreaks=120,
        ):
        """
        Writes sequences to file in phylip format, interleaving optional
        If outfile = 'stdout' the sequences are printed to screen, not written
        to disk
        If print_to_screen = True the sequences are printed to screen
        whether they are written to disk or not
        """

        maxlen = len(max(self.sequences, key=len))
        file_header = ' {0} {1}'.format(self.length, maxlen)
        s = [file_header]
        maxheader = len(max(self.headers, key=len))
        label_length = max(maxheader + 1, 10)
        if interleaved:
            seq_length = linebreaks - label_length
            num_lines = maxlen / seq_length
            if maxlen % seq_length:
                num_lines += 1

            for i in range(num_lines):
                for seq_header in self.headers:
                    if i == 0:
                        s.append('{0:<{1}} {2}'.format(seq_header,
                                 label_length,
                                 (self.mapping[seq_header])[i
                                 * seq_length:(i + 1) * seq_length]))
                    else:
                        s.append('{0} {1}'.format(' ' * label_length,
                                 (self.mapping[seq_header])[i
                                 * seq_length:(i + 1) * seq_length]))
                s.append('')
        else:
            lines = ['{0:<{1}} {2:-<{3}}'.format(x, label_length, y,
                     maxlen) for (x, y) in zip(self.headers,
                     self.sequences)]
            s.extend(lines)
            s.append('')
        s = '\n'.join(s)
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
            return outfile


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
        headers=[],
        sequences=[],
        dv=[],
        datatype=None,
        tree=None,
        ):

        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.mapping = {}
        self.datatype = datatype
        self.length = 0
        self.seqlength = 0
        self.is_aligned = False
        self.TCfiles = {}
        self.dv = dv
        if tree:
            if isinstance(tree, Tree):
                self.tree = tree
        else:
            self.tree = Tree()
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
        if items == []:
            return self
        sort_key = lambda item: tuple((int(num) if num else alpha)
                for (num, alpha) in re.findall(r'(\d+)|(\D+)', item[0]))
        items = sorted(items, key=sort_key)
        (h, s) = zip(*items)
        if in_place:
            self.headers = h
            self.sequences = s
            return self
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

    def hashname(self):
        H = hashlib.sha1()
        H.update(self.name)
        return H.hexdigest()

    def _write_temp_phylip(self, tmpdir='/tmp', use_hashname=False):
        if use_hashname:
            filename = self.hashname()
        else:
            filename = self.name
        self.write_phylip('{0}/{1}.phy'.format(tmpdir, filename))
        return filename

    def _write_temp_tc(
        self,
        tmpdir='/tmp',
        make_guide_tree=True,
        use_hashname=True,
        ):
        if use_hashname:
            filename = self.hashname()
        else:
            filename = self.name
        num_matrices = len(self.dv)
        if num_matrices == 0:
            print 'No distance-variance matrix available'
            return
        all_labels = self.headers
        len_all_labels = len(all_labels)

        # Temporary files

        dv_tmpfile = open('{0}/{1}_dv.txt'.format(tmpdir, filename), 'w'
                          )
        labels_tmpfile = open('{0}/{1}_labels.txt'.format(tmpdir,
                              filename), 'w')
        map_tmpfile = open('{0}/{1}_map.txt'.format(tmpdir, filename),
                           'w')
        if make_guide_tree:
            guidetree_tmpfile = open('{0}/{1}_tree.nwk'.format(tmpdir,
                    filename), 'w')

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

        if make_guide_tree:
            guidetree = self.get_bionj_tree(ncat=1, tmpdir=tmpdir)
            dpy_guidetree = dpy.Tree()
            dpy_guidetree.read_from_string(guidetree.newick, 'newick')
            dpy_guidetree.resolve_polytomies()
            newick_string = dpy_guidetree.as_newick_string() + ';\n'
            guidetree_tmpfile.write(newick_string)
            guidetree_tmpfile.flush()
            guidetree_tmpfile.close()

        # Check it all worked

        assert os.path.isfile('{0}/{1}_dv.txt'.format(tmpdir, filename))
        assert os.path.isfile('{0}/{1}_labels.txt'.format(tmpdir,
                              filename))
        assert os.path.isfile('{0}/{1}_map.txt'.format(tmpdir,
                              filename))
        if make_guide_tree:
            assert os.path.isfile('{0}/{1}_tree.nwk'.format(tmpdir,
                                  filename))
        return filename

    def get_phyml_tree(
        self,
        model=None,
        datatype=None,
        ncat=4,
        tmpdir='/tmp',
        overwrite=True,
        verbose=False,
        ):

        if not overwrite and self.tree.newick:
            print '{0}: Tree exists and overwrite set to false'.format(self.name)
            return self.tree
        self.tree = Tree()
        filename = self._write_temp_phylip(tmpdir=tmpdir, use_hashname=True)
        print 'Running phyml on ' + str(self.name) + '...'
        input_file = '{0}/{1}.phy'.format(tmpdir, filename)
        if not model and not datatype:  # quick-fix to allow specification of other
            if self.datatype == 'dna':  # models when calling phyml
                model = 'GTR'
                datatype = 'nt'
            elif self.datatype == 'protein':
                model = 'WAG'
                datatype = 'aa'
            else:
                print 'I don\'t know this datatype: {0}'.format(self.datatype)
                return
        t = self.tree.run_phyml(
            model,
            input_file,
            datatype,
            self.name,
            ncat=ncat,
            overwrite=overwrite,
            verbose=verbose,
            )
        os.remove('{0}/{1}.phy'.format(tmpdir, filename))
        return self.tree

    def get_bionj_tree(
        self,
        model=None,
        datatype=None,
        ncat=1,
        optimise='n',
        tmpdir='/tmp',
        overwrite=True,
        verbose=False,
        ):

        if not overwrite and self.tree.newick:
            print '{0}: Tree exists and overwrite set to false'.format(self.name)
            return self.tree
        self.Tree = Tree()
        filename = self._write_temp_phylip(tmpdir=tmpdir, use_hashname=True)
        print 'Running bionj on ' + str(self.name) + '...'
        input_file = '{0}/{1}.phy'.format(tmpdir, filename)

        if not model and not datatype:  # quick-fix to allow specification of other
            if self.datatype == 'dna':  # models when calling phyml
                model = 'GTR'
                datatype = 'nt'
            elif self.datatype == 'protein':
                model = 'WAG'
                datatype = 'aa'
            else:
                print 'I don\'t know this datatype: {0}'.format(self.datatype)
                return
        t = self.tree.run_bionj(
            model,
            input_file,
            datatype,
            ncat=ncat,
            name=self.name,
            optimise=optimise,
            overwrite=overwrite,
            verbose=verbose,
            )
        os.remove('{0}/{1}.phy'.format(tmpdir, filename))
        return self.tree

    def get_raxml_tree(self, tmpdir='/tmp', overwrite=True):
        if not overwrite and self.tree.newick:
            print '{0}: Tree exists and overwrite set to false'.format(self.name)
            return self.tree
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
        self.tree.run_raxml(model, input_file, self.name, tmpdir,
                            overwrite=overwrite)
        os.remove('{0}/{1}.phy'.format(tmpdir, self.name))
        if os.path.isfile('{0}/{1}.phy.reduced'.format(tmpdir,
                          self.name)):
            os.remove('{0}/{1}.phy.reduced'.format(tmpdir, self.name))
        return self.tree

    def get_guide_tree(self, tmpdir='/tmp', overwrite=True):
        filename = self._write_temp_phylip(tmpdir=tmpdir, use_hashname=True)
        input_file = '{0}/{1}.phy'.format(tmpdir, filename)
        if self.datatype == 'dna':
            model = 'GTRGAMMA'
        elif self.datatype == 'protein':
            model = 'PROTGAMMAWAG'
        else:
            print 'I don\'t know this datatype: {0}'.format(self.datatype)
            return
        t = Tree().run_raxml(model, input_file, self.name, tmpdir,
                             guide=True)
        if os.path.isfile('{0}/{1}.phy'.format(tmpdir, filename)):
            os.remove('{0}/{1}.phy'.format(tmpdir, filename))
        if os.path.isfile('{0}/{1}.phy.reduced'.format(tmpdir,
                          filename)):
            os.remove('{0}/{1}.phy.reduced'.format(tmpdir, filename))
        return t

    def get_TC_tree(self, tmpdir='/tmp', overwrite=True):
        if not overwrite and self.tree.newick:
            print '{0}: Tree exists and overwrite set to false'.format(self.name)
            return self.tree
        filename = self._write_temp_tc(tmpdir=tmpdir, use_hashname=True)
        print 'Running TreeCollection on ' + str(self.name) + '...'
        self.tree.run_treecollection(
            '{0}/{1}_dv.txt'.format(tmpdir, filename),
            '{0}/{1}_map.txt'.format(tmpdir, filename),
            '{0}/{1}_labels.txt'.format(tmpdir, filename),
            '{0}/{1}_tree.nwk'.format(tmpdir, filename),
            self.name,
            overwrite=overwrite,
            )
        os.remove('{0}/{1}_dv.txt'.format(tmpdir, filename))
        os.remove('{0}/{1}_map.txt'.format(tmpdir, filename))
        os.remove('{0}/{1}_labels.txt'.format(tmpdir, filename))
        os.remove('{0}/{1}_tree.nwk'.format(tmpdir, filename))
        return self.tree

    def get_dv_matrix(
        self,
        tmpdir='/tmp',
        helper='/Users/kgori/Projects/clustering_project/class_files/DV_wrapper.drw'
            ,
        overwrite=True,
        ):
        """
        Makes a call to the TC_wrapper.drw darwin helper script, which
        calculates a distance-variance matrix from the sequence alignments,
        and generates files needed by the treecollection binary
        """

        if not overwrite and self.dv:
            return self.dv[0]
        if self.name:
            fastafile = '{0}/{1}.fas'.format(tmpdir, self.name)
        else:
            fastafile = '{0}/fasta_tmp.fas'.format(tmpdir)
        if self.datatype == 'dna':
            datatype = 'DNA'
        else:
            datatype = 'AA'

        if not os.path.isfile(helper):
            print 'Can\'t find the darwin helper file at {0}'.format(helper)
            return 0

        self.write_fasta(fastafile)
        print 'Running darwin on {0}, datatype = {1}'.format(fastafile,
                datatype, helper)
        command = \
            'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\'; fpath := \'{2}/\'; ReadProgram(\'{3}\');" | darwin'.format(fastafile,
                datatype, tmpdir, helper)
        # print command
        process = Popen(command, shell=True, stdout=PIPE, stderr=PIPE)
        (stdout, stderr) = process.communicate()
        dv_string = \
            open('{0}/temp_distvar.txt'.format(tmpdir)).read().rstrip()
        labels = ' '.join(self.headers)
        os.remove(fastafile)
        os.remove('{0}/temp_distvar.txt'.format(tmpdir))
        return (dv_string, labels)

    def _pivot(self, lst):
        new_lst = zip(*lst)
        return [''.join(x) for x in new_lst]

    def shuffle(self):
        """
        Modifies in-place
        """
        columns = self._pivot(self.sequences)
        shf(columns)
        self.sequences = self._pivot(columns)
        self._update()

    def split_by_lengths(self, lengths, names=None):
        assert sum(lengths) == self.seqlength
        columns = self._pivot(self.sequences)
        newcols = []
        for l in lengths:
            newcols.append(columns[:l])
            columns = columns[l:]
        newrecs = []
        for col in newcols:
            newseqs = self._pivot(col)
            newrec = TCSeqRec(headers=self.headers,
                              sequences=newseqs, datatype=self.datatype)
            newrecs.append(newrec)
        if names:
            for i, newrec in enumerate(newrecs):            
                newrec.name = names[i]
        else:
            for i, newrec in enumerate(newrecs):            
                newrec.name = 'record_{0}'.format(i+1)
        return newrecs
