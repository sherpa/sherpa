#!/usr/bin/env python

#from distutils.core import setup, Extension
from numpy.distutils.core import setup, Extension

setup(name='group', version='4.5', ext_modules=[
    Extension('group',
              ['pygrplib.c'],
              ['../src'],
              library_dirs=['../src/.libs/'],
              libraries=['grp'],
              depends=['pygrplib.h']
             )]
    )

