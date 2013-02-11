#!/usr/bin/env python

from . import ExternalSoftware
from ..errors import FileError, filecheck_and_raise
from ..utils.file_utils import *

local_dir = path_to(__file__)


class DV_wrapper(ExternalSoftware):

    default_binary = 'darwin'
    DV_wrapper = locate_file('DV_wrapper.drw', directory=local_dir)

    def __str__(self):
        desc = 'Wrapper for darwin method DV_wrapper.drw'

    def call(self, filename, datatype):
        datatype = datatype or 'AA'
        cmd = \
            'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\'; fpath := \'{2}/\'; ReadProgram(\'{3}\');" | {4}'.format(filename,
                datatype, self.tmpdir, self.DV_wrapper, self.binary)

        return subprocess(cmd)

    def read_datatype(self, sequence_record):
        if sequence_record.datatype == 'protein':
            return 'AA'
        elif sequence_record.datatype == 'dna':
            return 'DNA'

    def read(self):
        output = '{0}/temp_distvar.txt'.format(self.tmpdir)
        filecheck_and_raise(output)
        self.add_tempfile(output)
        return open(output).read().rstrip()

    def run(self, sequence_record):
        filename = self.write(sequence_record)
        datatype = self.read_datatype(sequence_record)
        (stdout, stderr) = self.call(filename, datatype)
        dv_string = self.read()
        labels = ' '.join(sequence_record.headers)
        self.clean()
        return (dv_string, labels)

    def write(self, sequence_record):
        if sequence_record.name:
            filename = (sequence_record.name if len(sequence_record.name)
                        < 40 else sequence_record.hashname())
        else:
            filename = 'fasta_tmp'
        filename = '{0}/{1}.tmp'.format(self.tmpdir, filename)
        sequence_record.write_fasta(filename)
        self.add_tempfile(filename)
        return filename
