"""C++ wrapper for pyne library header."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector

cdef extern from "pyne.h" namespace "pyne":
    std_string PYNE_DATA
    std_string NUC_DATA_PATH

    double endftod(char *) except +
    void pyne_start() except +
