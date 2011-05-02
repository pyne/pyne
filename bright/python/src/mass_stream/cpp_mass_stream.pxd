"""C++ wrapper for isoname library."""
from libcpp.map cimport map
from libcpp.set cimport set

cimport std


cdef extern from "../../../cpp/MassStream.h":
    cdef cppclass MassStream:
        # Constuctors
        MassStream()
        MassStream(map[int, double], float, std.string) except +
        MassStream(char *, float, std.string) except +

        # Attributes
        map[int, double] comp
        double mass
        std.string name

        # Methods
        void norm_comp_dict() except +
        void load_from_hdf5(std.string, std.string, int) except +
        void load_from_text(char *) except +

        void print_ms() except +
        void normalize() except +
        map[int, double] mult_by_mass() except +
        double atomic_weight() except +

        # Substream Methods
        MassStream get_sub_stream(set[int], std.string) except +
        #MassStream get_sub_stream(set[std.string], std.string) except + # Redundant
        MassStream get_u(std.string) except +
        MassStream get_pu(std.string) except +
        MassStream get_lan(std.string) except +
        MassStream get_act(std.string) except +
        MassStream get_tru(std.string) except +
        MassStream get_ma(std.string) except +
        MassStream get_fp(std.string) except +

        # Operator Overloads
        MassStream operator+(double) except +
        MassStream operator+(MassStream) except +
        MassStream operator*(double) except +
        MassStream operator/(double) except +
