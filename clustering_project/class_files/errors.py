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


class FilenameError(FileError):

    """
    Raise when trying to make a file or directory
    with an invalid name - in Unix the only forbidden
    character is '/'
    """

    def __str__(self):
        return '\'{0}\' is an invalid filename'.format(self.value)


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


def directorycheck_and_make(directory):
    try:
        directorycheck_and_raise(directory)
    except DirectoryError, e:
        print e
        try:
            if '/' in directory:
                raise FilenameError(directory)
            print 'Creating \'{0}\''.format(directory)
            os.mkdir(directory)
        except FilenameError, e:
            print e


def directorycheck_and_make_recursive(directory):
    try:
        directorycheck_and_raise(directory)
    except DirectoryError, e:
        print e
        print 'Creating \'{0}\''.format(directory)
        os.makedirs(directory)
