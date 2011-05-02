# Cython imports
from libcpp.map cimport map as cpp_map
from cython import pointer


# Local imports
cimport std
cimport cpp_mass_stream

cdef class MassStream:
    cdef cpp_mass_stream.MassStream * ms_pointer


# Dictionary - Map Converters
ctypedef cpp_mass_stream.MassStream * msp
cdef cpp_map[std.string, msp] dict_to_map_str_msp(dict)
cdef dict map_to_dict_str_msp(cpp_map[std.string, msp])
