#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension

rast = Extension('_rast',
                 library_dirs=['.'],
                 libraries = ['rast'],
                 swig_opts = ["-c++"],
                 sources=['rast.i'])

setup (name = 'rast',
       version = '0.0',
       author      = "Thomas Breuel",
       description = """rast library bindings""",
       ext_modules = [rast],
       py_modules = ["rast"],
       )
