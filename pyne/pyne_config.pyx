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

PYNE_CORE_BASE_PATH = os.path.join(__path__[0], "core")

# TODO: Remove this when Scikit build support is enabled 
PYNE_CORE_BASE_PATH = os.path.join(__path__[0], '..', '..', '..', '..')

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


def get_paths(subdir, pattern="*", recursive=False):
    """
    Helper function to return paths that match a given pattern within a subdirectory.

    Args:
        subdir (str): The subdirectory within the 'core' directory.
        pattern (str): The pattern to match files or directories.
        recursive (bool): Whether to search recursively in subdirectories.

    Returns:
        list: A list of matched paths.
    """
    search_pattern = (
        os.path.join(PYNE_CORE_BASE_PATH, subdir, "**", pattern)
        if recursive
        else os.path.join(PYNE_CORE_BASE_PATH, subdir, pattern)
    )
    return glob.glob(search_pattern, recursive=recursive)


def get_include_path():
    """Return includes and include path for PyNE headers."""
    include = get_paths("include", "*", recursive=True)
    include_path = get_paths("include", "", recursive=False)
    return include, include_path


def get_core_libraries():
    """Return libraries and library paths for PyNE."""
    lib = [
        lib_file
        for lib in ["lib", "lib64"]
        for lib_file in get_paths(lib, "libpyne*", recursive=True)
    ]
    lib_path = [
        lib_file
        for lib in ["lib", "lib64"]
        for lib_file in get_paths(lib, "", recursive=False)
    ]
    return lib, lib_path


def get_extra_libraries():
    """Return the extra libraries installed by auditwheel or delocate."""
    libs_path = (
        os.path.join(__path__[0], ".dylibs")
        if sys.platform == "darwin"
        else os.path.normpath(os.path.join(__path__[0], "..", "pyne.libs"))
    )
    return (
        (glob.glob(os.path.join(libs_path, "*")), libs_path)
        if os.path.exists(libs_path)
        else ([], [])
    )


# Setup variables
include, include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()

# Export variables for easy access
__all__ = ["include", "include_path", "lib", "lib_path", "extra_lib", "extra_lib_path"]
