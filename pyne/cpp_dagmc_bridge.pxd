"""C++ wrapper for dagmc_bridge."""
from libcpp.string cimport string as std_string
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.utility cimport pair

cimport extra_types

cdef extern from "moab/Types.hpp" namespace "moab":

    ctypedef size_t EntityHandle
    ctypedef enum ErrorCode:
        pass

cdef extern from "dagmc_bridge.h" namespace "pyne":

    ctypedef double vec3[3]

    int dag_ent_handle_size() except +
    float dag_version() except +
    extra_types.uint32 dag_rev_version() except +
    const int * geom_id_list(int, int*) except +
    EntityHandle handle_from_id(int dimension, int id) except +
    int id_from_handle(EntityHandle eh) except +
    ErrorCode dag_load(const char* filename) except + 
    ErrorCode dag_pt_in_vol(EntityHandle vol, vec3 pt, int* result, vec3 dir,
                            const void* history) except +
