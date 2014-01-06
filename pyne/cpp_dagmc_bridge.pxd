"""C++ wrapper for dagmc_bridge."""
from libcpp.string cimport string as std_string
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.utility cimport pair

cimport extra_types

cdef extern from "dagmc_bridge.h" namespace "pyne":

    int dag_ent_handle_size(void) except +
    float dag_version(void) except +
    extra_types.uint32 dag_rev_version(void) except +
