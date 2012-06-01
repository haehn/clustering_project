#!/usr/bin/python
# -*- coding: utf-8 -*-


class SequenceRecord(object):

    """ 
    Class for reading sequence files in fasta or phylip formats
    Supports writing in fasta, phylip, nexus formats,
    sorting sequences by length and name,
    concatenating sequences when sequence names are a perfect match,
    iterating over records.
    """

    def get_fasta_file(self, fasta_file, name=None):
        """ FASTA format parser: turns fasta file into Alignment_record object """

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
        self._update()

    def get_phylip_file(self, phylip_file, name=None):
        """ PHYLIP format parser"""

        headers = []
        sequences = []
        openfile = open(phylip_file, 'r')
        info = openfile.readline().split()
        num_taxa = info[0]
        seq_length = info[1]

        while True:
            line = openfile.readline().rstrip()
            if not line:
                break
            line = line.split()
            headers.append(line[0])
            sequences.append(line[1])

        self.name = name
        (self.headers, self.sequences) = (headers, sequences)
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
        if infile:
            if file_format == 'fasta':
                self.get_fasta_file(infile, name)
            elif file_format == 'phylip':
                self.get_phylip_file(infile, name)
        self.index = -1
        self._update()

    def _update(self):
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
            output_string += \
                '''>{0}
{1}...({2})
'''.format(self.headers[i],
                    (self.sequences[i])[:50], len(self.sequences[i]))
        output_string += '{0} sequences in record'.format(len(self))
        return output_string

    def __len__(self):
        return self.length

    def __add__(self, other):

        try:
            assert set(self.headers) == set(other.headers)
        except AssertionError:
            print 'Sequence labels do not match between alignments'
            return self
        d = {}
        for k in self.mapping.keys():
            d[k] = self.mapping[k] + other.mapping[k]
        return SequenceRecord(headers=d.keys(), sequences=d.values())

    def __mul__(self, n):
        if not isinstance(n, int):
            print 'Truncating {0} to {1}'.format(n, int(n))
            n = int(n)
        new_seqs = []
        for s in self.sequences:
            new_seqs.append(s * n)
        return SequenceRecord(headers=self.headers, sequences=new_seqs)

    def sort_by_length(self):

        # Sort sequences by descending order of length
        # Uses zip as its own inverse [ zip(*zip(A,B)) == (A,B) ]

        (h, s) = zip(*sorted(zip(self.headers, self.sequences),
                     key=lambda item: len(item[1].replace('-', '')),
                     reverse=True))
        self.headers = h
        self.sequences = s
        return self

    def sort_by_name(self):
        (h, s) = zip(*sorted(zip(self.headers, self.sequences)))
        self.headers = h
        self.sequences = s
        return self

    def write_fasta(self, outfile='stdout', print_to_screen=False):
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
        file_header = '''#NEXUS

'''
        file_header += 'begin data;\n'
        file_header += \
            '    dimensions ntax={0} nchar={1};\n'.format(self.length,
                maxlen)
        file_header += \
            '    format datatype={0} interleave=no gap=-;\n'.format(sequence_type)
        file_header += '''    matrix

'''

        file_footer = '''    ;

end;
'''

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
