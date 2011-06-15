#!/usr/bin/env python
 
import os
from copy import deepcopy

from distutils.core import setup
from distutils.extension import Extension
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from Cython.Distutils import build_ext

import numpy as np

from setup_data import INFO


##########################################
### Set compiler options for extensions ###
###########################################
#pyt_dir = os.path.abspath(os.path.join('pyne'))
#cpp_dir = os.path.abspath(os.path.join('..', 'cpp'))
#dat_dir = os.path.abspath(os.path.join('..', 'data'))

pyt_dir = os.path.join('pyne')
cpp_dir = os.path.join('..', 'cpp')
dat_dir = os.path.join('..', 'data')


# Get numpy include dir
numpy_include = np.get_include()



#
# For stlconverters
# 
stlconv_ext = {'name': "pyne.stlconverters"}

stlconv_ext['sources'] = [os.path.join('pyne', 'stlconverters.pyx')]
stlconv_ext['include_dirs'] = [pyt_dir, cpp_dir, numpy_include]
stlconv_ext['language'] = "c++"

if os.name == 'posix':
#    stlconv_ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
#    stlconv_ext["undef_macros"] = ["NDEBUG"]
    pass
elif os.name == 'nt':
    stlconv_ext["extra_compile_args"] = ["/EHsc"]
    stlconv_ext["define_macros"] = [("_WIN32", None)]


#
# For nucname
# 
nucname_ext = {'name': 'pyne.nucname'}

nucname_ext['sources'] = [
    'pyne.cpp',
    'nucname.cpp',
    ]
nucname_ext['sources'] = [os.path.join(cpp_dir, s) for s in nucname_ext['sources']] + \
                         [os.path.join('pyne', 'nucname.pyx')]

nucname_ext['include_dirs'] = [pyt_dir, cpp_dir, numpy_include]
nucname_ext['language'] = "c++"


if os.name == 'posix':
    #isoname_ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
    #nucname_ext["undef_macros"] = ["NDEBUG"]
    pass
elif os.name == 'nt':
    nucname_ext["extra_compile_args"] = ["/EHsc"]
    nucname_ext["define_macros"] = [("_WIN32", None)]




##########################
### Setup Package Data ###
##########################
pack_dir = {
    'pyne': os.path.join('pyne'), 
    }
    
pack_data = {'pyne': []}


###################
### Call setup! ###
###################
setup(name="pyne",
    version = INFO['version'],
    description = 'Python for Nuclear Engineering',
    author = 'PyNE Development Team',
    author_email = 'scopatz@gmail.com',
    url = 'https://github.com/pyne',
    packages = ['pyne'],
    package_dir = pack_dir,
    cmdclass = {'build_ext': build_ext}, 
    ext_modules=[
        Extension(**stlconv_ext), 
        Extension(**nucname_ext), 
        ],
    )

