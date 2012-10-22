#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
import numpy as np

ext_modules = [Extension('evrot_extensions', ['evrot_extensions.c'], 
    include_dirs=[np.get_include()])]

setup(name='evrot_extensions',
    version='1.0',
    ext_modules=ext_modules,
    include_dirs=[np.get_include()]
)
