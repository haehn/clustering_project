#!/usr/bin/env python

import os
import sys


class FileError(Exception):

    """
    Reporting when file errors occur
    """

    def __init__(self, filename=''):
        """
        Store the file that gave the exception
        """

        self.value = filename

    def __str__(self):
        return 'Error opening file {0}'.format(self.value)


class DirectoryError(FileError):

    """
    Reporting when directory errors occur
    """

    def __str__(self):
        line1 = 'Error opening directory \'{0}\''.format(self.value)
        line2 = 'Directory doesn\'t exist'
        return '\n'.join((line1, line2))


def filecheck_and_raise(filename):
    if not os.path.isfile(filename):
        raise FileError(filename)


def filecheck_and_quit(filename):
    try:
        filecheck_and_raise(filename)
    except FileError, e:
        print e
        sys.exit()


def directorycheck_and_raise(directory):
    if not os.path.isdir(directory):
        raise DirectoryError(directory)


def directorycheck_and_quit(directory):
    try:
        if not os.path.isdir(directory):
            raise DirectoryError(directory)
    except DirectoryError, e:
        print e
        sys.exit()


def directorycheck_and_make(directory, verbose=True):
    try:
        directorycheck_and_raise(directory)
    except DirectoryError, e:
        if verbose:
            print e
            print 'Creating \'{0}\''.format(directory)
        os.makedirs(directory)
