import os
import sys
from ctypes import cdll

lib = os.path.join(os.path.split(__file__)[0], 'lib')
libtype = '.so'
if os.name == 'nt':
    libtype = '.dll'
elif sys.platform == 'darwin':
    libtype = '.dylib'
liblist = ['libpyne', 'libpyne_data', 'libpyne_material', 'libpyne_enrichemnt', 'libpyne_nucname', 'libpyne_rxname']
for item in liblist:
    cdll.LoadLibrary(os.path.join(lib, item + libtype))

from .pyne_config import *
