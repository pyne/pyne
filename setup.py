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


# Thanks to http://patorjk.com/software/taag/  
# and http://www.chris.com/ascii/index.php?art=creatures/dragons
# for ASCII art inspiriation

pyne_logo = """\

                                  /   \       
 _                        )      ((   ))     (                          
(@)                      /|\      ))_((     /|\                          
|-|                     / | \    (/\|/\)   / | \                      (@) 
| | -------------------/--|-voV---\`|'/--Vov-|--\---------------------|-|
|-|                         '^`   (o o)  '^`                          | |
| |                               `\Y/'                               |-|
|-|                                                                   | |
| |        /\             ___           __  __             /\         |-|
|-|       /^~\           / _ \_   _  /\ \ \/__\           /^~\        | |  
| |       /^~\          / /_)/ | | |/  \/ /_\             /^~\        |-|
|-|       /^~\         / ___/| |_| / /\  //__             /^~\        | | 
| |       ^||`         \/     \__, \_\ \/\__/             ^||`        |-|  
|-|        ||                |____/                        ||         | | 
| |       ====                                            ====        |-|
|-|                                                                   | |
| |                                                                   |-|
|-|___________________________________________________________________| |
(@)              l   /\ /         ( (       \ /\   l                `\|-|
                 l /   V           \ \       V   \ l                  (@)
                 l/                _) )_          \I                   
                                   `\ /'
                                     `  
"""


###########################################
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
    #ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}']
    ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}', '${ORIGIN}/.']
    #ext['runtime_library_dirs'] = ['${ORIGIN}/lib', '${ORIGIN}', '${ORIGIN}/.'] + \
    #                              [os.path.abspath(p) for p in ext['library_dirs']] + \
    #                              [os.path.abspath(p + '/pyne/lib') for p in sys.path] + \
    #                              [os.path.abspath(p + '/pyne') for p in sys.path] + \
    #                              [os.path.abspath(p) for p in sys.path]

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
exts.append(cpp_ext("pyne.lib.libpyne_nucname", ['nucname.cpp'], ['pyne']))

# data
exts.append(cpp_ext("pyne.lib.libpyne_data", ['data.cpp'], ['pyne', 'pyne_nucname'], True))

# material
exts.append(cpp_ext("pyne.lib.libpyne_material", ['material.cpp'], ['pyne', 'pyne_nucname', 'pyne_data'], True))


# Python extension modules
# STL converters
exts.append(cpp_ext("pyne.stlconverters", ['stlconverters.pyx']))

# pyne_config
exts.append(cpp_ext("pyne.pyne_config", ['pyne_config.pyx'], ['pyne']))

# nucname
exts.append(cpp_ext("pyne.nucname", ['nucname.pyx'], ['pyne', 'pyne_nucname']))

# material
exts.append(cpp_ext("pyne.data", ['data.pyx'], ['pyne', 'pyne_nucname', 'pyne_data'], True))

# material
exts.append(cpp_ext("pyne.material", ['material.pyx'], ['pyne', 'pyne_nucname', 'pyne_data', 'pyne_material'], True))



##########################
### Setup Package Data ###
##########################
packages = ['pyne', 'pyne.lib', 'pyne.dbgen']

pack_dir = {'pyne': 'pyne', 'pyne.dbgen': 'pyne/dbgen'}

pack_data = {'pyne': ['includes/*.h', 'includes/pyne/*.pxd'],
             'pyne.dbgen': ['*.html'],
            }

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
    print pyne_logo

    # clean includes dir and recopy files over
    if os.path.exists('pyne/includes'):
        remove_tree('pyne/includes')

    mkpath('pyne/includes')
    for header in glob.glob('cpp/*.h'):
        copy_file(header, 'pyne/includes')

    mkpath('pyne/includes/pyne')
    for header in glob.glob('pyne/*.pxd'):
        copy_file(header, 'pyne/includes/pyne')

    # call setup
    setup(name="pyne",
        version = INFO['version'],
        description = 'Python for Nuclear Engineering',
        author = 'PyNE Development Team',
        author_email = 'scopatz@gmail.com',
        url = 'http://pyne.github.com/pyne',
        packages = packages,
        package_dir = pack_dir,
        package_data = pack_data,
        cmdclass = {'build_ext': build_ext}, 
        ext_modules=ext_modules,
        scripts=scripts, 
        )

    # Clean includes after setup has run
    if os.path.exists('pyne/includes'):
        remove_tree('pyne/includes')
