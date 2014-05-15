

from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "endf.h" namespace "pyne::endf":
    cdef struct mt_451_struct:
        int nuc_id
        double awr
        int mat
        int mf
        int mt
        int lrp
        int lfi
        int nlib
        int nmod
        double elis
        int sta
        int lis
        int liso
        int nfor
        double awi
        double emax
        int lrel
        int nsub
        int nver
        double temp
        int ldrv
        vector[vector[int]] mt_list
    cdef struct endf_library_struct:
        mt_451_struct mt_451
    
    endf_library_struct library
    
    void read_endf(string) except +