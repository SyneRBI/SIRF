#!/usr/bin/env python

"""
setup.py file for SWIG SIRFReg
"""

from distutils.core import setup, Extension

SIRFReg_module = Extension('_SIRFReg', sources = ['SIRFReg.i'], swig_opts = ['-c++'])

setup (name        = 'SIRFReg',
       version     = '0.1',
       author      = "Richard Brown",
       description = """Python module of SIRFReg""",
       ext_modules = [SIRFReg_module],
       py_modules  = ["SIRFReg"]
       )