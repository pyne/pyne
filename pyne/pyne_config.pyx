"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport free

# local imports 
cimport std
cimport cpp_pyne

import os

nuc_data = os.path.join(os.path.split(__file__)[0], 'nuc_data.h5')


####################################
### pyne configuration namespace ###
####################################

# Expose the C-code start up routine
def pyne_start():
    # Specifiy the BRIGHT_DATA directory
    if "PYNE_DATA" not in os.environ:
        bd = os.path.split(__file__)
        os.environ['PYNE_DATA'] = os.path.join(*(bd[0], ''))

    # Specifiy the NUC_DATA_PATH 
    if "NUC_DATA_PATH" not in os.environ:
        os.environ['NUC_DATA_PATH'] = nuc_data

    # Call the C-version of pyne_start
    cpp_pyne.pyne_start()


# Run the appropriate start-up routines
pyne_start()


################################
### PyNE Configuration Class ###
################################
cdef class PyneConfig:
    """A PyNE configuration helper class."""

    property PYNE_DATA:
        def __get__(self):
            cdef std.string value = cpp_bright.PYNE_DATA
            return value.c_str()

        def __set__(self, char * value):
            cpp_bright.PYNE_DATA = std.string(value)


    property NUC_DATA_PATH:
        def __get__(self):
            cdef std.string value = cpp_bright.NUC_DATA_PATH
            return value.c_str()

        def __set__(self, char * value):
            cpp_bright.NUC_DATA_PATH = std.string(value)


        
# Make a singleton of the pyne config object
pyne_config = PyneConfig()


