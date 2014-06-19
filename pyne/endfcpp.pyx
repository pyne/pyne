from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map
from libcpp.vector cimport vector

cimport pyne.cpp_endf

cdef class _mt_base_struct:
    cdef pyne.cpp_endf.mt_base_struct *thisptr
    def __cinit__(self):
        pass
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:

cdef class _mt_451_struct:
    cdef pyne.cpp_endf.mt_451_struct *thisptr
    cdef pyne.cpp_endf.mt_451_struct ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:
    property nwd:
        def __get__(self): return self.thisptr.nwd
    property ldrv:
        def __get__(self): return self.thisptr.ldrv
    property sta:
        def __get__(self): return self.thisptr.sta
    property temp:
        def __get__(self): return self.thisptr.temp
    property emax:
        def __get__(self): return self.thisptr.emax
    property nlib:
        def __get__(self): return self.thisptr.nlib
    property lfi:
        def __get__(self): return self.thisptr.lfi
    property mt_list:
        def __get__(self): return self.thisptr.mt_list
    property nmod:
        def __get__(self): return self.thisptr.nmod
    property nfor:
        def __get__(self): return self.thisptr.nfor
    property lrp:
        def __get__(self): return self.thisptr.lrp
    property lis:
        def __get__(self): return self.thisptr.lis
    property awi:
        def __get__(self): return self.thisptr.awi
    property elis:
        def __get__(self): return self.thisptr.elis
    property nxc:
        def __get__(self): return self.thisptr.nxc
    property liso:
        def __get__(self): return self.thisptr.liso
    property nsub:
        def __get__(self): return self.thisptr.nsub
    property lrel:
        def __get__(self): return self.thisptr.lrel
    property nver:
        def __get__(self): return self.thisptr.nver

cdef class _mt_fpy_8_struct:
    cdef pyne.cpp_endf.mt_fpy_8_struct *thisptr
    cdef pyne.cpp_endf.mt_fpy_8_struct ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:
    property i:
        def __get__(self): return self.thisptr.i
    property yields:
        def __get__(self): return self.thisptr.yields
    property le:
        def __get__(self): return self.thisptr.le
    property e:
        def __get__(self): return self.thisptr.e

cdef class _mt_452_1_struct:
    cdef pyne.cpp_endf.mt_452_1_struct *thisptr
    cdef pyne.cpp_endf.mt_452_1_struct ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:
    property eint:
        def __get__(self): return self.thisptr.eint
    property nu_e:
        def __get__(self): return self.thisptr.nu_e
    property lnu:
        def __get__(self): return self.thisptr.lnu
    property intn:
        def __get__(self): return self.thisptr.intn
    property poly:
        def __get__(self): return self.thisptr.poly
    property nbt:
        def __get__(self): return self.thisptr.nbt

cdef class _mt_455_1_struct:
    cdef pyne.cpp_endf.mt_455_1_struct *thisptr
    cdef pyne.cpp_endf.mt_455_1_struct ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:
    property ldg:
        def __get__(self): return self.thisptr.ldg
    property nu_d:
        def __get__(self): return self.thisptr.nu_d
    property lnu:
        def __get__(self): return self.thisptr.lnu
    property intn:
        def __get__(self): return self.thisptr.intn
    property einti:
        def __get__(self): return self.thisptr.einti
    property ne:
        def __get__(self): return self.thisptr.ne
    property nbt:
        def __get__(self): return self.thisptr.nbt
    property alpha_arr:
        def __get__(self): return self.thisptr.alpha_arr
    property lambda_arr:
        def __get__(self): return self.thisptr.lambda_arr
    property eint:
        def __get__(self): return self.thisptr.eint
    property lambdas:
        def __get__(self): return self.thisptr.lambdas

cdef class _mt_456_1_struct:
    cdef pyne.cpp_endf.mt_456_1_struct *thisptr
    cdef pyne.cpp_endf.mt_456_1_struct ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:
    property eint:
        def __get__(self): return self.thisptr.eint
    property nu_e:
        def __get__(self): return self.thisptr.nu_e
    property lnu:
        def __get__(self): return self.thisptr.lnu
    property intn:
        def __get__(self): return self.thisptr.intn
    property nbt:
        def __get__(self): return self.thisptr.nbt
    property nu:
        def __get__(self): return self.thisptr.nu

cdef class _mt_458_1_struct:
    cdef pyne.cpp_endf.mt_458_1_struct *thisptr
    cdef pyne.cpp_endf.mt_458_1_struct ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:
    property degp:
        def __get__(self): return self.thisptr.degp
    property defr:
        def __get__(self): return self.thisptr.defr
    property egp:
        def __get__(self): return self.thisptr.egp
    property end:
        def __get__(self): return self.thisptr.end
    property efr:
        def __get__(self): return self.thisptr.efr
    property denp:
        def __get__(self): return self.thisptr.denp
    property det:
        def __get__(self): return self.thisptr.det
    property denu:
        def __get__(self): return self.thisptr.denu
    property eb:
        def __get__(self): return self.thisptr.eb
    property enp:
        def __get__(self): return self.thisptr.enp
    property degd:
        def __get__(self): return self.thisptr.degd
    property enu:
        def __get__(self): return self.thisptr.enu
    property et:
        def __get__(self): return self.thisptr.et
    property deb:
        def __get__(self): return self.thisptr.deb
    property egd:
        def __get__(self): return self.thisptr.egd
    property er:
        def __get__(self): return self.thisptr.er
    property der:
        def __get__(self): return self.thisptr.der
    property dend:
        def __get__(self): return self.thisptr.dend

cdef class _mt_460_1_struct:
    cdef pyne.cpp_endf.mt_460_1_struct *thisptr
    cdef pyne.cpp_endf.mt_460_1_struct ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        def __get__(self): return self.thisptr.mt
    property mf:
        def __get__(self): return self.thisptr.mf
    property awr:
        def __get__(self): return self.thisptr.awr
    property mat:
        def __get__(self): return self.thisptr.mat
    #property {}:
    property tint:
        def __get__(self): return self.thisptr.tint
    property elist:
        def __get__(self): return self.thisptr.elist
    property lo:
        def __get__(self): return self.thisptr.lo
    property intn:
        def __get__(self): return self.thisptr.intn
    property nbt:
        def __get__(self): return self.thisptr.nbt
    property t:
        def __get__(self): return self.thisptr.t
    property ng:
        def __get__(self): return self.thisptr.ng
    property lambdas:
        def __get__(self): return self.thisptr.lambdas


cdef class endf:
    cdef pyne.cpp_endf.library *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new pyne.cpp_endf.library()
    def __dealloc__(self):
        del self.thisptr
    def read_endf(self, filenm):
        """
        Read in an endf file and load the contents into memory

        Parameters
        ----------
        filename : str
            filename of an ENDF 6 formatted file
        """
        self.thisptr.read_endf(filenm)
    def get_content_list(self):
        """
        Get a list of the loaded contents

        Returns
        -------
        contents : vector[mat,mf,mt]
            Returns a vector of vectors containing a list of the current
            ENDF data loaded into memory
        """
        return self.thisptr.get_content_list()
    def get(self, mat, mf, mt):
        """
        Get a single record from the ENDF data in memory

        Parameters
        ----------
        mat : int
            material of interest in ENDF 6 format
        mf : int
            ENDF 6 file number
        mt : int
            ENDF 6 reaction number

        Returns
        -------
        mt_data : mt_struct
           Returns a python wrapper object of the struct of interest
        """
        if mf == 1:
            if mt == 451:
                mt_451 = _mt_451_struct()
                mt_451.ob = self.thisptr.get[pyne.cpp_endf.mt_451_struct](mat, mf, mt)
                return mt_451
            elif mt == 452:
                mt_452 = _mt_452_1_struct()
                mt_452.ob = self.thisptr.get[pyne.cpp_endf.mt_452_1_struct](mat, mf, mt)
                return mt_452
            elif mt == 455:
                mt_455 = _mt_455_1_struct()
                mt_455.ob = self.thisptr.get[pyne.cpp_endf.mt_455_1_struct](mat, mf, mt)
                return mt_455
            elif mt == 456:
                mt_456 = _mt_456_1_struct()
                mt_456.ob = self.thisptr.get[pyne.cpp_endf.mt_456_1_struct](mat, mf, mt)
                return mt_456
            elif mt == 458:
                mt_458 = _mt_458_1_struct()
                mt_458.ob = self.thisptr.get[pyne.cpp_endf.mt_458_1_struct](mat, mf, mt)
                return mt_458
            elif mt == 460:
                mt_460 = _mt_460_1_struct()
                mt_460.ob = self.thisptr.get[pyne.cpp_endf.mt_460_1_struct](mat, mf, mt)
                return mt_460
        elif mf == 8:
            if mt == 454 or mt == 459:
                mt_fpy = _mt_fpy_8_struct()
                mt_fpy.ob = (self.thisptr.get[pyne.cpp_endf.mt_fpy_8_struct](mat, mf, mt))
                return mt_fpy

    def get_mt(self, mf, mt):
        """
        Gets a list of all data in the library with a given mf and mt

        Parameters
        ----------
        mf : int
            ENDF 6 file number
        mt : int
            ENDF 6 reaction number

        Returns
        -------
        mt_list : list of mt_struct
            returns a list of structs with the desired mf and mt number

        """
        retlist = []
        cdef vector[pyne.cpp_endf.mt_451_struct] d_451_list
        cdef vector[pyne.cpp_endf.mt_452_1_struct] d_452_1_list
        cdef vector[pyne.cpp_endf.mt_455_1_struct] d_455_1_list
        cdef vector[pyne.cpp_endf.mt_456_1_struct] d_456_1_list
        cdef vector[pyne.cpp_endf.mt_458_1_struct] d_458_1_list
        cdef vector[pyne.cpp_endf.mt_460_1_struct] d_460_1_list
        cdef vector[pyne.cpp_endf.mt_fpy_8_struct] d_fpy_8_list

        if mf == 1:
            if mt == 451:
                d_451_list = self.thisptr.getl[pyne.cpp_endf.mt_451_struct](mf, mt)
                for i in range(d_451_list.size()):
                    mt_451 = _mt_451_struct()
                    mt_451.ob = d_451_list[i]
                    retlist.append(mt_451)
                return retlist
            elif mt == 452:
                d_452_1_list = self.thisptr.getl[pyne.cpp_endf.mt_452_1_struct](mf, mt)
                for i in range(d_452_1_list.size()):
                    mt_452 = _mt_452_1_struct()
                    mt_452.ob = d_452_1_list[i]
                    retlist.append(mt_452)
                return retlist
            elif mt == 455:
                d_455_1_list = self.thisptr.getl[pyne.cpp_endf.mt_455_1_struct](mf, mt)
                for i in range(d_455_1_list.size()):
                    mt_455 = _mt_455_1_struct()
                    mt_455.ob = d_455_1_list[i]
                    retlist.append(mt_455)
                return retlist
            elif mt == 456:
                d_456_1_list = self.thisptr.getl[pyne.cpp_endf.mt_456_1_struct](mf, mt)
                for i in range(d_456_1_list.size()):
                    mt_456 = _mt_456_1_struct()
                    mt_456.ob = d_456_1_list[i]
                    retlist.append(mt_456)
                return retlist
            elif mt == 458:
                d_458_1_list = self.thisptr.getl[pyne.cpp_endf.mt_458_1_struct](mf, mt)
                for i in range(d_458_1_list.size()):
                    mt_458 = _mt_458_1_struct()
                    mt_458.ob = d_458_1_list[i]
                    retlist.append(mt_458)
                return retlist
            elif mt == 460:
                d_460_1_list = self.thisptr.getl[pyne.cpp_endf.mt_460_1_struct](mf, mt)
                for i in range(d_460_1_list.size()):
                    mt_460 = _mt_460_1_struct()
                    mt_460.ob = d_460_1_list[i]
                    retlist.append(mt_460)
                return retlist
        elif mf == 8:
            if mt == 454 or mt == 459:
                d_fpy_8_list = self.thisptr.getl[pyne.cpp_endf.mt_fpy_8_struct](mf, mt)
                for i in range(d_fpy_8_list.size()):
                    mt_fpy = _mt_fpy_8_struct()
                    mt_fpy.ob = d_fpy_8_list[i]
                    retlist.append(mt_fpy)
                return retlist
