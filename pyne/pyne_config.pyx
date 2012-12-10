"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

# local imports 
include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
cimport cpp_pyne

import os
import json
import pyne.__init__

prefix = os.path.split(pyne.__init__.__file__)[0]

lib = os.path.join(prefix, 'lib')
includes = os.path.join(prefix, 'include')
nuc_data = os.path.join(prefix, 'nuc_data.h5')


####################################
### pyne configuration namespace ###
####################################

# Expose the C-code start up routine
def pyne_start():
    # Specifiy the BRIGHT_DATA directory
    if "PYNE_DATA" not in os.environ:
        os.environ['PYNE_DATA'] = prefix

    # Specifiy the NUC_DATA_PATH 
    if "NUC_DATA_PATH" not in os.environ:
        os.environ['NUC_DATA_PATH'] = nuc_data

    # load cached metadata
    with open(os.path.join(prefix, "metadata.json"), 'r') as f:
        md = json.load(f)

    # set HDF5 dir on path
    if os.getenv('HDF5_DIR'):
        md['HDF5_DIR'] = os.environ['HDF5_DIR']
    libdll = 'dll' if os.name == 'nt' else 'lib'
    ldpath = 'PATH' if os.name == 'nt' else 'LD_LIBRARY_PATH'
    sepcha = ';' if os.name == 'nt' else ':'
    #if isinstance(md['HDF5_DIR'], basestring) and 0 < len(md['HDF5_DIR']):
    #    os.environ[ldpath] += sepcha + os.path.join(md['HDF5_DIR'], libdll)
    
    # Call the C-version of pyne_start
    cpp_pyne.pyne_start()


# Run the appropriate start-up routines
pyne_start()


################################
### PyNE Configuration Class ###
################################
cdef class PyneConf:
    """A PyNE configuration helper class."""

    property PYNE_DATA:
        def __get__(self):
            cdef std_string value = cpp_pyne.PYNE_DATA
            return <char *> value.c_str()

        def __set__(self, char * value):
            cpp_pyne.PYNE_DATA = std_string(value)


    property NUC_DATA_PATH:
        def __get__(self):
            cdef std_string value = cpp_pyne.NUC_DATA_PATH
            return <char *> value.c_str()

        def __set__(self, char * value):
            cpp_pyne.NUC_DATA_PATH = std_string(value)


        
# Make a singleton of the pyne config object
pyne_conf = PyneConf()

# hacks for not communicating environment (windows issue)
if pyne_conf.PYNE_DATA == "<NOT_FOUND>":
    pyne_conf.PYNE_DATA = os.environ['PYNE_DATA']
if pyne_conf.NUC_DATA_PATH == "<NOT_FOUND>":
    pyne_conf.NUC_DATA_PATH = os.environ['NUC_DATA_PATH']


