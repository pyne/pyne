#ifndef PYNE_SKQ36P4BFNE3VI6VHVADCDT4VQ
#define PYNE_SKQ36P4BFNE3VI6VHVADCDT4VQ

/* Types.hpp may be validly included from C code */
#include <moab/Types.hpp>
#include <DagMC.hpp>

#ifdef __cplusplus
using moab::ErrorCode;
using moab::EntityHandle;
using moab::DagMC;

namespace pyne {

static DagMC* DAG = new DagMC();
  
extern "C" {
#endif
/* Notice to future maintainers:
 * If this file is ever used as a C header from a C compilation unit, 
 * there will probably need to be an #elseif here that defines EntityHandle
 * or includes an appropriate MOAB header.  I don't think Types.hpp guarantees
 * that EntityHandle will be defined in C.
 */

typedef double vec3[3];

float dag_version(void);

unsigned dag_rev_version(void);

int dag_ent_handle_size(void);

const int* geom_id_list(int dimension, int* number_of_items);

EntityHandle handle_from_id(int dimension, int id);
int id_from_handle(EntityHandle eh);

ErrorCode dag_load(const char* filename);

void* dag_alloc_ray_history(void);

void dag_dealloc_ray_history(void* history);

ErrorCode dag_ray_fire(EntityHandle vol, vec3 ray_start, vec3 ray_dir, 
                        EntityHandle* next_surf, double* next_surf_dist,
                        void* history, double distance_limit);

ErrorCode dag_ray_follow(EntityHandle firstvol, vec3 ray_start, vec3 ray_dir,
                          double distance_limit, int* num_intersections,
                          EntityHandle** surfs, double** distances, 
                          EntityHandle** volumes, void* data_buffers);

void dag_dealloc_ray_buffer(void* data_buffers);

ErrorCode dag_pt_in_vol(EntityHandle vol, vec3 pt, int* result, vec3 dir, 
                         const void* history);

ErrorCode dag_next_vol(EntityHandle surface, EntityHandle volume, 
                        EntityHandle* next_vol);

int vol_is_graveyard(EntityHandle vol);
/* int surf_is_spec_refl(EntityHandle surf); */
/* int surf_is_white_refl(EntityHandle surf); */
int vol_is_implicit_complement(EntityHandle vol);

ErrorCode get_volume_metadata(EntityHandle vol, int* material, double* density, double* importance);

ErrorCode get_volume_boundary(EntityHandle vol, vec3 minPt, vec3 maxPt);

#ifdef __cplusplus
} // namespace pyne
} // extern "C"
#endif

#endif /* PYNE_SKQ36P4BFNE3VI6VHVADCDT4VQ */
