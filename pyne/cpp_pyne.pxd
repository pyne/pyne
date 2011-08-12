"""C++ wrapper for pyne library header."""
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.vector cimport vector

cimport std


cdef extern from "../cpp/pyne.h" namespace "pyne":
    std.string PYNE_DATA
    std.string NUC_DATA_PATH

    void pyne_start() except +
