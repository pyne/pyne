"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../../../cpp/MassStream.h":
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

        void print_ms()
        void normalize()
        map[int, double] mult_by_mass()
        double atomic_weight()

        # Substream Methods
        MassStream get_sub_stream(set[int], std.string)
        #MassStream get_sub_stream(set[std.string], std.string) # Redundant
        MassStream get_u(std.string)
        MassStream get_pu(std.string)
        MassStream get_lan(std.string)
        MassStream get_act(std.string)
        MassStream get_tru(std.string)
        MassStream get_ma(std.string)
        MassStream get_fp(std.string)

        # Operator Overloads
        MassStream operator+(double)
        MassStream operator+(MassStream)
        MassStream operator*(double)
        MassStream operator/(double)
