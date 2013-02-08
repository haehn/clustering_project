#!/usr/bin/env python


class ExternalSoftware(object):

    def __init__(self, binary, tmpdir):
        pass

    def __str__(self):
        pass

    def add_flag(self, flag, value):
        self.flags[flag] = value

    def call(self):
        pass

    def clean(self):
        pass

    def read(self):
        pass

    def run(self):
        pass

    def writetmp(self):
        pass


class Phyml(ExternalSoftware):

    pass
