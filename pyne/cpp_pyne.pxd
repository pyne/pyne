"""C++ wrapper for pyne library header."""
from libcpp.map cimport map
from libcpp.set cimport set

include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
    from libcpp.vector cimport vector
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
    from pyne._includes.libcpp.vector cimport vector

cdef extern from "pyne.h" namespace "pyne":
    std_string PYNE_DATA
    std_string NUC_DATA_PATH

    double endftod(char *) except +
    void pyne_start() except +
