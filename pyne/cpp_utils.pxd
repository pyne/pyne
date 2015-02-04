"""C++ wrapper for utils header."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector

cdef extern from "utils.h" namespace "pyne":
    std_string PYNE_DATA
    std_string NUC_DATA_PATH
    bint USE_WARNINGS 
    bint toggle_warnings() except +

    double endftod(char *) except +
    void use_fast_endftod() except +
    void pyne_start() except +
