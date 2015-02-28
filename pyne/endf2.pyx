from libcpp.vector cimport vector
cimport pyne.cpp_endf

cdef class _mt_base:
    cdef pyne.cpp_endf.mt_base *thisptr
    cdef pyne.cpp_endf.mt_base ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:

cdef class _mt451:
    cdef pyne.cpp_endf.mt451 *thisptr
    cdef pyne.cpp_endf.mt451 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property nwd:
        """ number of descriptove records"""
        def __get__(self): return self.thisptr.nwd
    property ldrv:
        """ derived material flag"""
        def __get__(self): return self.thisptr.ldrv
    property sta:
        """ target stability"""
        def __get__(self): return self.thisptr.sta
    property temp:
        """ temperature in kelvin"""
        def __get__(self): return self.thisptr.temp
    property emax:
        """ upper energy limit"""
        def __get__(self): return self.thisptr.emax
    property nlib:
        """ library identifier"""
        def __get__(self): return self.thisptr.nlib
    property lfi:
        """ does material fission"""
        def __get__(self): return self.thisptr.lfi
    property mt_list:
        """ list of [mf,mt,lines,mod]"""
        def __get__(self): return self.thisptr.mt_list
    property nmod:
        """ material modification"""
        def __get__(self): return self.thisptr.nmod
    property nfor:
        """ library format"""
        def __get__(self): return self.thisptr.nfor
    property lrp:
        """ resonance parameters given"""
        def __get__(self): return self.thisptr.lrp
    property lis:
        """ state number of the nucleas"""
        def __get__(self): return self.thisptr.lis
    property awi:
        """ projectile mass relative to neutron"""
        def __get__(self): return self.thisptr.awi
    property elis:
        """ excitation energy"""
        def __get__(self): return self.thisptr.elis
    property nxc:
        """ number of mt record in this file"""
        def __get__(self): return self.thisptr.nxc
    property liso:
        """ isomeric state number"""
        def __get__(self): return self.thisptr.liso
    property nsub:
        """ sub library number"""
        def __get__(self): return self.thisptr.nsub
    property lrel:
        """ library release number"""
        def __get__(self): return self.thisptr.lrel
    property nver:
        """ library version number"""
        def __get__(self): return self.thisptr.nver

cdef class _mtfpy_mf8:
    cdef pyne.cpp_endf.mtfpy_mf8 *thisptr
    cdef pyne.cpp_endf.mtfpy_mf8 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property i:
        """ interpolation to be used between E[i-1] and E[i]"""
        def __get__(self): return self.thisptr.i
    property yields:
        """ yield data [zafp,fps,yi,dyi]"""
        def __get__(self): return self.thisptr.yields
    property le:
        """ number of energy  dependent yields given"""
        def __get__(self): return self.thisptr.le
    property e:
        """ list of energies"""
        def __get__(self): return self.thisptr.e

cdef class _mt452_mf1:
    cdef pyne.cpp_endf.mt452_mf1 *thisptr
    cdef pyne.cpp_endf.mt452_mf1 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property eint:
        """ Energy of the incident neutron"""
        def __get__(self): return self.thisptr.eint
    property nu_e:
        """ Neutrons per fission at the given energy"""
        def __get__(self): return self.thisptr.nu_e
    property lnu:
        """ type of data in section"""
        def __get__(self): return self.thisptr.lnu
    property intn:
        """ list of interpolations to be used"""
        def __get__(self): return self.thisptr.intn
    property poly:
        """ polynomial describing neutrons per fission(E)"""
        def __get__(self): return self.thisptr.poly
    property nbt:
        """ list of interpolation segments"""
        def __get__(self): return self.thisptr.nbt

cdef class _mt455_mf1:
    cdef pyne.cpp_endf.mt455_mf1 *thisptr
    cdef pyne.cpp_endf.mt455_mf1 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property ldg:
        """ energy dependence of decay constants"""
        def __get__(self): return self.thisptr.ldg
    property nu_d:
        """ average neutrons per delayed fission event"""
        def __get__(self): return self.thisptr.nu_d
    property lnu:
        """ data representation type"""
        def __get__(self): return self.thisptr.lnu
    property intn:
        """ list of interpolations to be used"""
        def __get__(self): return self.thisptr.intn
    property einti:
        """ energy interpolation scheme"""
        def __get__(self): return self.thisptr.einti
    property ne:
        """ number of energies"""
        def __get__(self): return self.thisptr.ne
    property nbt:
        """ list of interpolation segments"""
        def __get__(self): return self.thisptr.nbt
    property alpha_arr:
        """"""
        def __get__(self): return self.thisptr.alpha_arr
    property lambda_arr:
        """"""
        def __get__(self): return self.thisptr.lambda_arr
    property eint:
        """ ith energy"""
        def __get__(self): return self.thisptr.eint
    property lambdas:
        """ decay constant of ith precursor"""
        def __get__(self): return self.thisptr.lambdas

cdef class _mt456_mf1:
    cdef pyne.cpp_endf.mt456_mf1 *thisptr
    cdef pyne.cpp_endf.mt456_mf1 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property eint:
        """ ith energy"""
        def __get__(self): return self.thisptr.eint
    property nu_e:
        """ average neutrons per prompt fission event"""
        def __get__(self): return self.thisptr.nu_e
    property lnu:
        """"""
        def __get__(self): return self.thisptr.lnu
    property intn:
        """ list of interpolations to be used"""
        def __get__(self): return self.thisptr.intn
    property nbt:
        """ list of interpolation segments"""
        def __get__(self): return self.thisptr.nbt
    property nu:
        """ polynomial describing neutrons per fission(E)"""
        def __get__(self): return self.thisptr.nu

cdef class _mt458_mf1:
    cdef pyne.cpp_endf.mt458_mf1 *thisptr
    cdef pyne.cpp_endf.mt458_mf1 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property efr:
        """ kinetic energy from fission products"""
        def __get__(self): return self.thisptr.efr
    property egp:
        """ total energy released by prompt gamma rays"""
        def __get__(self): return self.thisptr.egp
    property eb:
        """ total energy released by delayed betas"""
        def __get__(self): return self.thisptr.eb
    property pen:
        """ kinetic energy of prompt fission neutrons"""
        def __get__(self): return self.thisptr.pen
    property enu:
        """ energy carried away by neutrinos"""
        def __get__(self): return self.thisptr.enu
    property den:
        """ kinetic energy of delayed fission neutrons"""
        def __get__(self): return self.thisptr.den
    property et:
        """ total energy release per fission"""
        def __get__(self): return self.thisptr.et
    property er:
        """ total energy release not lost to neutrinos"""
        def __get__(self): return self.thisptr.er
    property egd:
        """ total energy released by delayed gamma rays"""
        def __get__(self): return self.thisptr.egd

cdef class _mt460_mf1:
    cdef pyne.cpp_endf.mt460_mf1 *thisptr
    cdef pyne.cpp_endf.mt460_mf1 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property tint:
        """ time of ith multiplicity"""
        def __get__(self): return self.thisptr.tint
    property elist:
        """ energy of the ith photon in eV"""
        def __get__(self): return self.thisptr.elist
    property lo:
        """ representation type: 1 if discrete, 2 if continuous"""
        def __get__(self): return self.thisptr.lo
    property intn:
        """ list^2 of interpolations"""
        def __get__(self): return self.thisptr.intn
    property nbt:
        """ list^2 interpolation segments"""
        def __get__(self): return self.thisptr.nbt
    property t:
        """ Time dependence of ith multiplicity"""
        def __get__(self): return self.thisptr.t
    property ng:
        """ number"""
        def __get__(self): return self.thisptr.ng
    property lambdas:
        """ decay constant for the ith precursor"""
        def __get__(self): return self.thisptr.lambdas

cdef class _mt457_mf8:
    cdef pyne.cpp_endf.mt457_mf8 *thisptr
    cdef pyne.cpp_endf.mt457_mf8 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property ricl:
        """ l-shell internal conversion coefficient"""
        def __get__(self): return self.thisptr.ricl
    property lcon:
        """ Continuum spectrum flag"""
        def __get__(self): return self.thisptr.lcon
    property rtyp:
        """ decay mode of the nuclide"""
        def __get__(self): return self.thisptr.rtyp
    property ner:
        """ number of tabulated discrete energies"""
        def __get__(self): return self.thisptr.ner
    property lis:
        """ state of the original nuclide"""
        def __get__(self): return self.thisptr.lis
    property nd:
        """ number of decay modes given """
        def __get__(self): return self.thisptr.nd
    property ris:
        """ internal pair formation coefficient"""
        def __get__(self): return self.thisptr.ris
    property styp:
        """ Decay radiation"""
        def __get__(self): return self.thisptr.styp
    property ricc:
        """ total internal conversion coefficient"""
        def __get__(self): return self.thisptr.ricc
    property type:
        """ type of beta/ec transition"""
        def __get__(self): return self.thisptr.type
    property fc:
        """ continuum spectrum normalization"""
        def __get__(self): return self.thisptr.fc
    property fd:
        """ discrete spectrum normalization factor"""
        def __get__(self): return self.thisptr.fd
    property rick:
        """ k-shell internal conversion coefficient"""
        def __get__(self): return self.thisptr.rick
    property liso:
        """ isomeric state of the original nuclide"""
        def __get__(self): return self.thisptr.liso
    property nsp:
        """ total number of radiation types for which information is given"""
        def __get__(self): return self.thisptr.nsp
    property er:
        """ discrete energy of radiation"""
        def __get__(self): return self.thisptr.er
    property nst:
        """ nucleus stability flag: 0 = radioactive"""
        def __get__(self): return self.thisptr.nst
    property ri:
        """ intensity of discrete radiation"""
        def __get__(self): return self.thisptr.ri
    property erel:
        """ energy released"""
        def __get__(self): return self.thisptr.erel
    property eav:
        """ average decay energy of radiation"""
        def __get__(self): return self.thisptr.eav

cdef class _mf3:
    cdef pyne.cpp_endf.mf3 *thisptr
    cdef pyne.cpp_endf.mf3 ob
    def __cinit__(self):
        self.thisptr = &self.ob
    def __dealloc__(self):
        pass
    property nuc_id:
        """ ENDF style nuc id ZZZAAA"""
        def __get__(self): return self.thisptr.nuc_id
    property mt:
        """ ENDF reaction designation"""
        def __get__(self): return self.thisptr.mt
    property mf:
        """ ENDF file number"""
        def __get__(self): return self.thisptr.mf
    property awr:
        """ mass relative to neutron [double]"""
        def __get__(self): return self.thisptr.awr
    property mat:
        """ ENDF material"""
        def __get__(self): return self.thisptr.mat
    #property {}:
    property energy:
        """ ith energy"""
        def __get__(self): return self.thisptr.energy
    property intn:
        """ list of interpolations to be used"""
        def __get__(self): return self.thisptr.intn
    property sigma:
        """ ith cross section value"""
        def __get__(self): return self.thisptr.sigma
    property nbt:
        """ list of interpolation segments"""
        def __get__(self): return self.thisptr.nbt


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
        return self.thisptr.content_list
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
        mtdata : mt
           Returns a python wrapper object of the struct of interest
        """
        if mf == 1:
            if mt == 451:
                mt451 = _mt451()
                mt451.ob = self.thisptr.get[pyne.cpp_endf.mt451](mat, mf, mt)
                return mt451
            elif mt == 452:
                mt452 = _mt452_mf1()
                mt452.ob = self.thisptr.get[pyne.cpp_endf.mt452_mf1](mat, mf, mt)
                return mt452
            elif mt == 455:
                mt455 = _mt455_mf1()
                mt455.ob = self.thisptr.get[pyne.cpp_endf.mt455_mf1](mat, mf, mt)
                return mt455
            elif mt == 456:
                mt456 = _mt456_mf1()
                mt456.ob = self.thisptr.get[pyne.cpp_endf.mt456_mf1](mat, mf, mt)
                return mt456
            elif mt == 458:
                mt458 = _mt458_mf1()
                mt458.ob = self.thisptr.get[pyne.cpp_endf.mt458_mf1](mat, mf, mt)
                return mt458
            elif mt == 460:
                mt460 = _mt460_mf1()
                mt460.ob = self.thisptr.get[pyne.cpp_endf.mt460_mf1](mat, mf, mt)
                return mt460
        elif mf == 8:
            if mt == 454 or mt == 459:
                mtfpy = _mtfpy_mf8()
                mtfpy.ob = (self.thisptr.get[pyne.cpp_endf.mtfpy_mf8](mat, mf, mt))
                return mtfpy
            elif mt == 457:
                mt457 = _mt457_mf8()
                mt457.ob = (self.thisptr.get[pyne.cpp_endf.mt457_mf8](mat, mf, mt))
                return mt457
        elif mf == 3:
            mf3 = _mf3()
            mf3.ob = (self.thisptr.get[pyne.cpp_endf.mf3](mat, mf, mt))
            return mf3

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
        mt_list : list of mt
            returns a list of structs with the desired mf and mt number

        """
        retlist = []
        cdef vector[pyne.cpp_endf.mt451] d_451_list
        cdef vector[pyne.cpp_endf.mt452_mf1] d_452_1_list
        cdef vector[pyne.cpp_endf.mt455_mf1] d_455_1_list
        cdef vector[pyne.cpp_endf.mt456_mf1] d_456_1_list
        cdef vector[pyne.cpp_endf.mt458_mf1] d_458_1_list
        cdef vector[pyne.cpp_endf.mt460_mf1] d_460_1_list
        cdef vector[pyne.cpp_endf.mtfpy_mf8] d_fpy_8_list
        cdef vector[pyne.cpp_endf.mt457_mf8] d_457_8_list
        cdef vector[pyne.cpp_endf.mf3] d_mf3_list

        if mf == 1:
            if mt == 451:
                d_451_list = self.thisptr.getl[pyne.cpp_endf.mt451](mf, mt)
                for i in range(d_451_list.size()):
                    mt451 = _mt451()
                    mt451.ob = d_451_list[i]
                    retlist.append(mt451)
                return retlist
            elif mt == 452:
                d_452_1_list = self.thisptr.getl[pyne.cpp_endf.mt452_mf1](mf, mt)
                for i in range(d_452_1_list.size()):
                    mt452 = _mt452_mf1()
                    mt452.ob = d_452_1_list[i]
                    retlist.append(mt452)
                return retlist
            elif mt == 455:
                d_455_1_list = self.thisptr.getl[pyne.cpp_endf.mt455_mf1](mf, mt)
                for i in range(d_455_1_list.size()):
                    mt455 = _mt455_mf1()
                    mt455.ob = d_455_1_list[i]
                    retlist.append(mt455)
                return retlist
            elif mt == 456:
                d_456_1_list = self.thisptr.getl[pyne.cpp_endf.mt456_mf1](mf, mt)
                for i in range(d_456_1_list.size()):
                    mt456 = _mt456_mf1()
                    mt456.ob = d_456_1_list[i]
                    retlist.append(mt456)
                return retlist
            elif mt == 458:
                d_458_1_list = self.thisptr.getl[pyne.cpp_endf.mt458_mf1](mf, mt)
                for i in range(d_458_1_list.size()):
                    mt458 = _mt458_mf1()
                    mt458.ob = d_458_1_list[i]
                    retlist.append(mt458)
                return retlist
            elif mt == 460:
                d_460_1_list = self.thisptr.getl[pyne.cpp_endf.mt460_mf1](mf, mt)
                for i in range(d_460_1_list.size()):
                    mt460 = _mt460_mf1()
                    mt460.ob = d_460_1_list[i]
                    retlist.append(mt460)
                return retlist
        elif mf == 8:
            if mt == 454 or mt == 459:
                d_fpy_8_list = self.thisptr.getl[pyne.cpp_endf.mtfpy_mf8](mf, mt)
                for i in range(d_fpy_8_list.size()):
                    mtfpy = _mtfpy_mf8()
                    mtfpy.ob = d_fpy_8_list[i]
                    retlist.append(mtfpy)
                return retlist
            elif mt == 457:
                d_457_8_list = self.thisptr.getl[pyne.cpp_endf.mt457_mf8](mf, mt)
                for i in range(d_457_8_list.size()):
                    mt457 = _mt457_mf8()
                    mt457.ob = d_457_8_list[i]
                    retlist.append(mt457)
                return retlist
        elif mf == 3:
            d_mf3_list = self.thisptr.getl[pyne.cpp_endf.mf3](mf, mt)
            for i in range(d_mf3_list.size()):
                mf3 = _mf3()
                mf3.ob = d_mf3_list[i]
                retlist.append(mf3)
            return retlist
