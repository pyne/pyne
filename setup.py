#!/usr/bin/env python
 
import os
import glob
from copy import deepcopy

from distutils.core import setup
from distutils.extension import Extension
from distutils.util import get_platform
from distutils.file_util import copy_file, move_file
from distutils.dir_util import mkpath, remove_tree
from distutils.sysconfig import get_python_version
from Cython.Distutils import build_ext

import numpy as np
import sys

INFO = {
    'version': '0.1',
    }



##########################################
### Set compiler options for extensions ###
###########################################
pyt_dir = os.path.join('pyne')
cpp_dir = os.path.join('cpp')
dat_dir = os.path.join('data')


# Get numpy include dir
numpy_include = np.get_include()

# HDF5 stuff
posix_hdf5_libs = ["z", "m", "hdf5", "hdf5_hl", "hdf5_cpp", "hdf5_hl_cpp",]
nt_hdf5_libs = ["/DEFAULTLIB:szip.lib", "/DEFAULTLIB:zlib1.lib", "/DEFAULTLIB:hdf5dll.lib",
                "/DEFAULTLIB:hdf5_hldll.lib", "/DEFAULTLIB:hdf5_cppdll.lib", "/DEFAULTLIB:hdf5_hl_cppdll.lib", ]
nt_hdf5_extra_compile_args = ["/EHsc"]
nt_hdf5_macros = [("_WIN32", None), ("_HDF5USEDLL_", None), ("HDF5CPP_USEDLL", None), ]


def cpp_ext(name, sources, libs=None, use_hdf5=False):
    """Helper function for setting up extension dictionary.

    Parameters
    ----------
    name : str
        Module name
    sources : list of str
        Files to compile
    libs : list of str
        Additional files to link against
    use_hdf5 : bool
        Link against hdf5?
    """
    ext = {'name': name}

    ext['sources'] = [os.path.join(cpp_dir, s) for s in sources if s.endswith('cpp')] + \
                     [os.path.join(pyt_dir, s) for s in sources if s.endswith('pyx')] + \
                     [s for s in sources if not any([s.endswith(suf) for suf in ['cpp', 'pyx']])]

    ext["libraries"] = []
    ext['include_dirs'] = [pyt_dir, cpp_dir, numpy_include]
    ext['language'] = "c++"

    # may need to be more general
    ext['library_dirs'] = ['build/lib/pyne/lib',
                           'build/lib.{0}-{1}/pyne/lib'.format(get_platform(), get_python_version()),
                           ]
    # perfectly general, thanks to dynamic runtime linking of $ORIGIN
    ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}']


    if os.name == 'posix':
        #ext["extra_compile_args"] = ["-Wno-strict-prototypes"]
        ext["undef_macros"] = ["NDEBUG"]
        if use_hdf5:
            ext["libraries"] += posix_hdf5_libs
        if libs is not None:
            ext["libraries"] += libs
    elif os.name == 'nt':
        ext["extra_compile_args"] = ["/EHsc"]
        ext["define_macros"] = [("_WIN32", None)]

        if use_hdf5:
            ext["libraries"] += nt_hdf5_libs
            ext["extra_compile_args"] += nt_hdf5_extra_compile_args
            ext["define_macros"] += nt_hdf5_macros

        if libs is not None:
            ext["libraries"] += libs

    return ext


#
# For extensions
# 
exts = []

# Pure C/C++ share libraries
# pyne lib
exts.append(cpp_ext("pyne.lib.libpyne", ['pyne.cpp']))

# nucname
exts.append(cpp_ext("pyne.lib.libpyne_nucname", ['nucname.cpp'], ['pyne'], True))

# material
exts.append(cpp_ext("pyne.lib.libpyne_material", ['material.cpp'], ['pyne', 'pyne_nucname'], True))


# Python extension modules
# STL converters
exts.append(cpp_ext("pyne.stlconverters", ['stlconverters.pyx']))

# pyne_config
exts.append(cpp_ext("pyne.pyne_config", ['pyne_config.pyx'], ['pyne']))

# nucname
exts.append(cpp_ext("pyne.nucname", ['nucname.pyx'], ['pyne', 'pyne_nucname'], True))

# material
exts.append(cpp_ext("pyne.material", ['material.pyx'], ['pyne', 'pyne_nucname', 'pyne_material'], True))



##########################
### Setup Package Data ###
##########################
pack_dir = {
    'pyne': os.path.join('pyne'), 
    }
    
pack_data = {'pyne': []}

ext_modules=[Extension(**ext) for ext in exts] 

# Compiler directives
compiler_directives = {'embedsignature': False}
for e in ext_modules:
    e.pyrex_directives = compiler_directives


# Utility scripts
scripts=['scripts/nuc_data_make']

###################
### Call setup! ###
###################
if __name__ == "__main__":
    setup(name="pyne",
        version = INFO['version'],
        description = 'Python for Nuclear Engineering',
        author = 'PyNE Development Team',
        author_email = 'scopatz@gmail.com',
        url = 'http://pyne.github.com/pyne',
        packages = ['pyne', 'pyne.lib', 'pyne.dbgen'],
        package_dir = pack_dir,
        cmdclass = {'build_ext': build_ext}, 
        ext_modules=ext_modules,
        scripts=scripts, 
        )

