"""Python wrapper for isoname library."""
from __future__ import unicode_literals
import os

# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free
from libcpp.string cimport string as std_string

# local imports
from pyne import __path__
cimport cpp_utils

nuc_data = os.path.join(__path__[0], 'nuc_data.h5')

####################################
### pyne configuration namespace ###
####################################

# Expose the C-code start up routine
def pyne_start():
    if "PYNE_DATA" not in os.environ:
        os.environ['PYNE_DATA'] = __path__[0]

    # Specifiy the NUC_DATA_PATH
    if "NUC_DATA_PATH" not in os.environ:
        os.environ['NUC_DATA_PATH'] = nuc_data

    libdll = 'dll' if os.name == 'nt' else 'lib'
    ldpath = 'PATH' if os.name == 'nt' else 'LD_LIBRARY_PATH'
    sepcha = ';' if os.name == 'nt' else ':'

    # Call the C-version of pyne_start
    cpp_utils.pyne_start()

# Run the appropriate start-up routines
pyne_start()

################################
### PyNE Configuration Class ###
################################
cdef class PyneConf:
    """A PyNE configuration helper class."""

    property PYNE_DATA:
        def __get__(self):
            cdef std_string value = cpp_utils.PYNE_DATA
            return <char *> value.c_str()

        def __set__(self, char * value):
            cpp_utils.PYNE_DATA = std_string(value)


    property NUC_DATA_PATH:
        def __get__(self):
            cdef std_string value = cpp_utils.NUC_DATA_PATH
            return <char *> value.c_str()

        def __set__(self, char * value):
            cpp_utils.NUC_DATA_PATH = std_string(value)


# Make a singleton of the pyne config object
pyne_conf = PyneConf()

# hacks for not communicating environment (windows issue)
if pyne_conf.PYNE_DATA == "<NOT_FOUND>":
    pyne_conf.PYNE_DATA = os.environ['PYNE_DATA']
if pyne_conf.NUC_DATA_PATH == "<NOT_FOUND>":
    pyne_conf.NUC_DATA_PATH = os.environ['NUC_DATA_PATH']
