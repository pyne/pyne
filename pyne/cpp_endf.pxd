from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.utility cimport pair
from libcpp.map cimport map


cdef extern from "endf_mt.h" namespace "pyne::endf":
    cdef cppclass mt_base_struct:
        int nuc_id
        int mt
        int mf
        double awr
        int mat
        # virtual ~mt_base_struct()
    cdef cppclass mt_451_struct(mt_base_struct):
        int nwd
        int ldrv
        int sta
        double temp
        double emax
        int nlib
        int lfi
        vector[vector[int]] mt_list
        int nmod
        int nfor
        int lrp
        int lis
        double awi
        double elis
        int nxc
        int liso
        int nsub
        int lrel
        int nver

    cdef cppclass mt_fpy_8_struct(mt_base_struct):
        vector[int] i
        vector[vector[vector[double]]] yields
        int le
        vector[double] e

    cdef cppclass mt_452_1_struct(mt_base_struct):
        vector[double] eint
        vector[double] nu_e
        int lnu
        vector[int] intn
        vector[double] poly
        vector[int] nbt

    cdef cppclass mt_455_1_struct(mt_base_struct):
        int ldg
        vector[double] nu_d
        int lnu
        vector[int] intn
        vector[int] einti
        vector[int] ne
        vector[int] nbt
        vector[vector[double]] alpha_arr
        vector[vector[double]] lambda_arr
        vector[double] eint
        vector[double] lambdas

    cdef cppclass mt_456_1_struct(mt_base_struct):
        vector[double] eint
        vector[double] nu_e
        int lnu
        vector[int] intn
        vector[int] nbt
        vector[double] nu

    cdef cppclass mt_458_1_struct(mt_base_struct):
        vector[double] degp
        vector[double] defr
        vector[double] egp
        vector[double] end
        vector[double] efr
        vector[double] denp
        vector[double] det
        vector[double] denu
        vector[double] eb
        vector[double] enp
        vector[double] degd
        vector[double] enu
        vector[double] et
        vector[double] deb
        vector[double] egd
        vector[double] er
        vector[double] der
        vector[double] dend

    cdef cppclass mt_460_1_struct(mt_base_struct):
        vector[vector[double]] tint
        vector[double] elist
        int lo
        vector[vector[int]] intn
        vector[vector[int]] nbt
        vector[vector[double]] t
        int ng
        vector[double] lambdas


cdef extern from "endf.h" namespace "pyne::endf":
    cdef cppclass endf_id_struct:
        int mat
        int mf
        int mt

    cdef cppclass library:
        map[endf_id_struct, mt_base_struct*] contents
        vector[vector[int]] get_content_list() except +
        void read_endf(string) except +
        T get[T](int,int,int) except +
        vector[T] getl[T](int, int) except +
