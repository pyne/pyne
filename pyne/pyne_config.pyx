"""Python wrapper for isoname library."""
from __future__ import unicode_literals
import os
import sys
import glob

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

def get_path(subdir, pattern="*"):
    """Helper function to return paths that match a given pattern within a subdirectory."""
    path = os.path.join(__path__[0], subdir)
    return glob.glob(os.path.join(path, pattern)) if os.path.exists(path) else []

def get_include_path():
    """Return the include directory path for PyNE headers."""
    return os.path.join(__path__[0], '..', '..', '..', '..', "include")

def get_core_libraries():
    """Return library paths and library directory paths."""
    #lib_paths = [os.path.join(sys.prefix, subdir) for subdir in ["lib", "lib64"]]
    # TODO: Temporary fix for old pyne
    lib_paths = [os.path.join(__path__[0], '..', '..', '..', '..', "lib")]
    libs = [lib for subdir in lib_paths for lib in get_path(subdir)]
    return libs, [path for path in lib_paths if os.path.exists(path)][0]

def get_extra_libraries():
    """List all the extra libraries of PyNE."""
    libs_path = os.path.join(__path__[0], ".dylibs") if sys.platform == "darwin" else os.path.join(__path__[0], "..", "pyne.libs")
    return (glob.glob(os.path.join(libs_path, "*")), libs_path) if os.path.exists(libs_path) else ([], [])

# Setup variables
include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()
nuc_data = os.path.join(__path__[0], 'nuc_data.h5')

# Export variables for easy access
__all__ = ["include_path", "lib", "lib_path", "extra_lib", "extra_lib_path", "nuc_data"]

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
