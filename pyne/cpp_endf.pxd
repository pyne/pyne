from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.string cimport string
from libcpp.map cimport map

cdef extern from "endf_mt.h" namespace "pyne::endf":
    cdef cppclass mt_base:
        int nuc_id # ENDF style nuc id ZZZAAA
        int mt # ENDF reaction designation
        int mf # ENDF file number
        double awr # mass relative to neutron [double]
        int mat # ENDF material
        # virtual ~mt_base() #

    cdef cppclass mt451(mt_base):
        int nwd # number of descriptove records
        int ldrv # derived material flag
        int sta # target stability
        double temp # temperature in kelvin
        double emax # upper energy limit
        int nlib # library identifier
        int lfi # does material fission
        vector[vector[int]] mt_list # list of [mf,mt,lines,mod]
        int nmod # material modification
        int nfor # library format
        int lrp # resonance parameters given
        int lis # state number of the nucleas
        double awi # projectile mass relative to neutron
        double elis # excitation energy
        int nxc # number of mt record in this file
        int liso # isomeric state number
        int nsub # sub library number
        int lrel # library release number
        int nver # library version number

    cdef cppclass mtfpy_mf8(mt_base):
        vector[int] i # interpolation to be used between E[i-1] and E[i]
        vector[vector[vector[double]]] yields # yield data [zafp,fps,yi,dyi]
        int le # number of energy  dependent yields given
        vector[double] e # list of energies

    cdef cppclass mt452_mf1(mt_base):
        vector[double] eint # Energy of the incident neutron
        vector[double] nu_e # Neutrons per fission at the given energy
        int lnu # type of data in section
        vector[int] intn # list of interpolations to be used
        vector[double] poly # polynomial describing neutrons per fission(E)
        vector[int] nbt # list of interpolation segments

    cdef cppclass mt455_mf1(mt_base):
        int ldg # energy dependence of decay constants
        vector[double] nu_d # average neutrons per delayed fission event
        int lnu # data representation type
        vector[int] intn # list of interpolations to be used
        vector[int] einti # energy interpolation scheme
        vector[int] ne # number of energies
        vector[int] nbt # list of interpolation segments
        vector[vector[double]] alpha_arr #
        vector[vector[double]] lambda_arr #
        vector[double] eint # ith energy
        vector[double] lambdas # decay constant of ith precursor

    cdef cppclass mt456_mf1(mt_base):
        vector[double] eint # ith energy
        vector[double] nu_e # average neutrons per prompt fission event
        int lnu #
        vector[int] intn # list of interpolations to be used
        vector[int] nbt # list of interpolation segments
        vector[double] nu # polynomial describing neutrons per fission(E)

    cdef cppclass mt458_mf1(mt_base):
        vector[pair[double,double]] efr # kinetic energy from fission products
        vector[pair[double,double]] egp # total energy released by prompt gamma rays
        vector[pair[double,double]] eb # total energy released by delayed betas
        vector[pair[double,double]] pen # kinetic energy of prompt fission neutrons
        vector[pair[double,double]] enu # energy carried away by neutrinos
        vector[pair[double,double]] den # kinetic energy of delayed fission neutrons
        vector[pair[double,double]] et # total energy release per fission
        vector[pair[double,double]] er # total energy release not lost to neutrinos
        vector[pair[double,double]] egd # total energy released by delayed gamma rays

    cdef cppclass mt460_mf1(mt_base):
        vector[vector[double]] tint # time of ith multiplicity
        vector[double] elist # energy of the ith photon in eV
        int lo # representation type: 1 if discrete, 2 if continuous
        vector[vector[int]] intn # list^2 of interpolations
        vector[vector[int]] nbt # list^2 interpolation segments
        vector[vector[double]] t # Time dependence of ith multiplicity
        int ng # number
        vector[double] lambdas # decay constant for the ith precursor

    cdef cppclass mt457_mf8(mt_base):
        vector[pair[double,double]] ricl # l-shell internal conversion coefficient
        vector[int] lcon # Continuum spectrum flag
        vector[double] rtyp # decay mode of the nuclide
        vector[int] ner # number of tabulated discrete energies
        int lis # state of the original nuclide
        vector[int] nd # number of decay modes given 
        vector[pair[double,double]] ris # internal pair formation coefficient
        vector[double] styp # Decay radiation
        vector[pair[double,double]] ricc # total internal conversion coefficient
        vector[double] type # type of beta/ec transition
        vector[pair[double,double]] fc # continuum spectrum normalization
        vector[pair[double,double]] fd # discrete spectrum normalization factor
        vector[pair[double,double]] rick # k-shell internal conversion coefficient
        int liso # isomeric state of the original nuclide
        int nsp # total number of radiation types for which information is given
        vector[pair[double,double]] er # discrete energy of radiation
        int nst # nucleus stability flag: 0 = radioactive
        vector[pair[double,double]] ri # intensity of discrete radiation
        vector[pair[double,double]] erel # energy released
        vector[pair[double,double]] eav # average decay energy of radiation

    cdef cppclass mf3(mt_base):
        vector[double] energy # ith energy
        vector[int] intn # list of interpolations to be used
        vector[double] sigma # ith cross section value
        vector[int] nbt # list of interpolation segments

cdef extern from "endf.h" namespace "pyne::endf":
    cdef cppclass endf_id:
        int mat
        int mf
        int mt

    cdef cppclass library:
        map[endf_id, mt_base*] contents
        vector[vector[int]] content_list
        void read_endf(string) except +
        T get[T](int,int,int) except +
        vector[T] getl[T](int, int) except +