"""C++ wrapper for dagmc_bridge."""
from libcpp.utility cimport pair
from libcpp.map cimport map
from libcpp.set cimport set
from libcpp.string cimport string as std_string

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
    void* dag_alloc_ray_history() except +
    void dag_dealloc_ray_history(void* history) except +
    void dag_dealloc_ray_buffer(void* data_buffers) except +
    ErrorCode dag_ray_fire(EntityHandle vol, vec3 ray_start, vec3 ray_dir,
                           EntityHandle* next_surf, double* next_surf_dist,
                           void* history, double distance_limit) except +
    ErrorCode dag_ray_follow(EntityHandle firstvol, vec3 ray_start, vec3 ray_dir,
                             double distance_limit, int* num_intersections,
                             EntityHandle** surfs, double** distances,
                             EntityHandle** volumes, void* data_buffers) except +
    ErrorCode dag_next_vol(EntityHandle surface, EntityHandle volume,
                           EntityHandle* next_vol) except +
    int vol_is_graveyard(EntityHandle vol) except +
    int vol_is_implicit_complement(EntityHandle vol) except +
    ErrorCode get_volume_metadata(EntityHandle vol, int* material, double* density, 
                                  double* importance) except +
    ErrorCode get_volume_boundary(EntityHandle vol, vec3 minPt, vec3 maxPt) except +
