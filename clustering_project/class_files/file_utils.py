#!/usr/bin/env python

import os
from subprocess import Popen, PIPE
from errors import FileError


def can_locate(filename):
    return (os.path.isfile(filename) if filename else False)


def can_open(directory):
    return (os.path.isdir(directory) if directory else False)


def path_to(filename):
    return os.path.dirname(filename)


def locate_by_env(filename, path=None):
    path = os.getenv(path) or os.getenv('PATH', os.defpath)
    for directory in path.split(os.pathsep):
        if verify(filename, directory):
            print directory
            return os.path.abspath(directory)
        f = locate_by_dir(filename, directory)
        if f:
            return f


def locate_by_dir(filename, directory=None):
    f = os.path.join(directory, filename)
    return (f if can_locate(f) else None)


def verify(filename, path):
    return can_locate(path) and os.path.basename(path) == filename


def join_path(*elements):
    return os.path.join(*elements)


def locate_file(filename, env_var='', directory=''):
    f = locate_by_env(filename, env_var) or locate_by_dir(filename,
            directory)
    return (f if can_locate(f) else None)


def syscall(cmd):
    return os.system(cmd)


def subprocess(cmd):
    process = Popen(cmd, shell=True, stdout=PIPE,
                               stderr=PIPE)
    return process.communicate()


def delete(filename):
    if can_locate(filename):
        return os.remove(filename)
    else:
        raise FileError(filename)
