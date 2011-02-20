"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../MassStream.h":
    cdef cppclass MassStream:
        # Constuctors
        MassStream()
        MassStream(map[int, double], float, std.string)
        MassStream(char *, float, std.string)

        # Attributes
        map[int, double] comp
        double mass
        std.string name

        # Methods
        void norm_comp_dict()
        void load_from_hdf5(std.string, std.string, int) except +
        void load_from_text(char *) except +

        void Print()
        void Normalize()
        map[int, double] multByMass()
        double atomic_weight()

        # Substream Methods
        MassStream getSubStream(set[int], std.string)
        #MassStream getSubStream(set[std.string], std.string) # Redundant
        MassStream getU(std.string)
        MassStream getPU(std.string)
        MassStream getLAN(std.string)
        MassStream getACT(std.string)
        MassStream getTRU(std.string)
        MassStream getMA(std.string)
        MassStream getFP(std.string)

        # Operator Overloads
        MassStream operator+(double)
        MassStream operator+(MassStream)
        MassStream operator*(double)
        MassStream operator/(double)
