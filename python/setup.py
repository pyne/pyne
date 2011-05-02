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

# Get numpy include dir
numpy_include = np.get_include()


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
#    ext_modules=[
#        Extension(**pyne_ext), 
#        ],
    )

