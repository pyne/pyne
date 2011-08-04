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

# HDF5 stuff
posix_hdf5_libs = ["z", "m", "hdf5", "hdf5_hl", "hdf5_cpp", "hdf5_hl_cpp",]
nt_hdf5_libs = ["/DEFAULTLIB:szip.lib", "/DEFAULTLIB:zlib1.lib", "/DEFAULTLIB:hdf5dll.lib",
                "/DEFAULTLIB:hdf5_hldll.lib", "/DEFAULTLIB:hdf5_cppdll.lib", "/DEFAULTLIB:hdf5_hl_cppdll.lib", ]
nt_hdf5_extra_compile_args = ["/EHsc"]
nt_hdf5_macros = [("_WIN32", None), ("_HDF5USEDLL_", None), ("HDF5CPP_USEDLL", None), ]


def cpp_ext(name, sources, use_hdf5=False):
    """He;lper function for setting up extension dictionary."""
    ext = {'name': name}

    ext['sources'] = [os.path.join(cpp_dir, s) for s in sources if s.endswith('cpp')] + \
                     [s for s in sources if not s.endswith('cpp')]

    ext['include_dirs'] = [pyt_dir, cpp_dir, numpy_include]
    ext['language'] = "c++"

    if os.name == 'posix':
        #ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
        #ext["undef_macros"] = ["NDEBUG"]
        if use_hdf5:
            ext["libraries"] = posix_hdf5_libs
    elif os.name == 'nt':
        ext["extra_compile_args"] = ["/EHsc"]
        ext["define_macros"] = [("_WIN32", None)]

        if use_hdf5:
            ext["libraries"] = nt_hdf5_libs
            ext["extra_compile_args"] += nt_hdf5_extra_compile_args
            ext["define_macros"] += nt_hdf5_macros

    return ext


#
# For extensions
# 
stlconv_ext = cpp_ext("pyne.stlconverters", [os.path.join('pyne', 'stlconverters.pyx')])

nucname_ext = cpp_ext("pyne.nucname", ['pyne.cpp',
                                       'nucname.cpp',
                                       os.path.join('pyne', 'nucname.pyx'),])

material_ext = cpp_ext("pyne.material", ['pyne.cpp',
                                         'nucname.cpp',
                                         'material.cpp',
                                         os.path.join('pyne', 'material.pyx'),], True)



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
        Extension(**material_ext), 
        ],
    )

