# Cython imports
from libcpp.map cimport map
from libcpp.vector cimport vector


cdef extern from "transmuters.h" namespace "pyne::transmuters":

    map[int, double] cram(vector[double], const map[int, double]) except +
    map[int, double] cram(vector[double], const map[int, double], const int) except +

