#ifndef PYNE_IS_AMALGAMATED
#include "dagmc_bridge.h"
#endif

#include <DagMC.hpp>
#include <moab/CartVect.hpp>

using moab::CartVect;

#include <vector>
#include <map>

using moab::DagMC;
using moab::EntityHandle;

#define CHECKERR(err) \
    if((err) != moab::MB_SUCCESS) return err;

namespace pyne {

float dag_version(void) {
    return DagMC::version();
}

unsigned dag_rev_version(void) {
    return DagMC::interface_revision();
}

int dag_ent_handle_size(void) {
    return sizeof(EntityHandle);
}

static std::vector<int> surfList;
static std::vector<int> volList;

const int* geom_id_list(int dimension, int* number_of_items) {
    switch(dimension) {
    case 2:
        *number_of_items = surfList.size();
        return &(surfList.front());
    case 3:
        *number_of_items = volList.size();
        return &(volList.front());
    default:
        *number_of_items = 0;
        return NULL;
    }
}

EntityHandle handle_from_id(int dimension, int id) {
    return DAG->entity_by_id(dimension, id);
}

int id_from_handle(EntityHandle eh) {
    return DAG->get_entity_id(eh);
}

ErrorCode dag_load(const char* filename){
    ErrorCode err;

    err = DAG->load_file(filename);
    CHECKERR(err);
    err = DAG->init_OBBTree();
    CHECKERR(err);

    std::vector<std::string> metadata_keys;
    metadata_keys.push_back("imp");
    metadata_keys.push_back("mat");
    metadata_keys.push_back("rho");
    metadata_keys.push_back("graveyard");

    std::map<std::string, std::string> metadata_synonyms;
    metadata_synonyms["rest.of.world"] = "graveyard";
    metadata_synonyms["outside.world"] = "graveyard";

    err = DAG->parse_properties(metadata_keys, metadata_synonyms);
    CHECKERR(err);

    int num_surfs = DAG->num_entities(2);
    surfList.reserve(num_surfs);
    for(int i = 1; i <= num_surfs; ++i) {
        surfList.push_back(DAG->id_by_index(2, i));
    }

    int num_vols = DAG->num_entities(3);
    volList.reserve(num_vols);
    for(int i = 1; i <= num_vols; ++i) {
        volList.push_back(DAG->id_by_index(3, i));
    }

    return err;
}


void* dag_alloc_ray_history(void) {
    return new DagMC::RayHistory();
}

void dag_dealloc_ray_history(void* r) {
    delete (static_cast<DagMC::RayHistory*>(r));
}

ErrorCode dag_ray_fire(EntityHandle vol, vec3 ray_start, vec3 ray_dir, 
                        EntityHandle* next_surf_ent, double* next_surf_dist,
                        void* history, double distance_limit) {
    ErrorCode err;

    DagMC*  dag = DAG;

    err = dag->ray_fire(vol, ray_start, ray_dir, *next_surf_ent, *next_surf_dist, 
                         static_cast<DagMC::RayHistory*>(history), distance_limit);
    CHECKERR(err);

    return err;
}

class ray_buffers {

    public:
    DagMC::RayHistory history;
    std::vector<EntityHandle> surfs;
    std::vector<double> dists;
    std::vector<EntityHandle> vols;

};

ErrorCode dag_ray_follow(EntityHandle firstvol, vec3 ray_start, vec3 ray_dir,
                          double distance_limit, int* num_intersections,
                          EntityHandle** surfs, double** distances, EntityHandle** volumes,
                          void* data_buffers){

    ray_buffers* buf = new ray_buffers;
    ErrorCode err;
    DagMC* dag = DAG;

    EntityHandle vol = firstvol;
    double dlimit = distance_limit;
    CartVect ray_point(ray_start);
    EntityHandle next_surf;
    double next_surf_dist;

    CartVect uvw(ray_dir);

    // iterate over the ray until no more intersections are available
    while(vol) {
        err = dag->ray_fire(vol, ray_point.array(), ray_dir, 
                             next_surf, next_surf_dist, &(buf->history), dlimit);
        CHECKERR(err);

        if(next_surf) {
            ray_point += uvw * next_surf_dist;
            buf->surfs.push_back(next_surf);
            buf->dists.push_back(next_surf_dist);
            err = dag->next_vol(next_surf, vol, vol);
            CHECKERR(err);
            buf->vols.push_back(vol);
            if(dlimit != 0){
                dlimit -= next_surf_dist;
            }
        }
        else vol = 0;
    }

    // assign to the output variables
    *num_intersections = buf->surfs.size();
    *surfs = &(buf->surfs[0]);
    *distances = &(buf->dists[0]);
    *volumes = &(buf->vols[0]);
    data_buffers = buf;

    return err;
}

void dag_dealloc_ray_buffer(void* data_buffers) {
    ray_buffers* b = static_cast<ray_buffers*>(data_buffers);
    delete b;
}

ErrorCode dag_pt_in_vol(EntityHandle vol, vec3 pt, int* result, vec3 dir, const void* history) {
    
    ErrorCode err;

    DagMC* dag = DAG;
    
    err = dag->point_in_volume(vol, pt, *result, dir, static_cast<const DagMC::RayHistory*>(history));

    return err;
}

ErrorCode dag_next_vol(EntityHandle surface, EntityHandle volume, EntityHandle* next_vol) {

    ErrorCode err;
    DagMC* dag = DAG;

    err = dag->next_vol(surface, volume, *next_vol);

    return err;
}

int vol_is_graveyard(EntityHandle vol) {
    return DAG->has_prop(vol, "graveyard");
}

/* int surf_is_spec_refl(EntityHandle surf); */
/* int surf_is_white_refl(EntityHandle surf); */

int vol_is_implicit_complement(EntityHandle vol){
    return DAG->is_implicit_complement(vol);
}

ErrorCode get_volume_metadata(EntityHandle vol, int* material, double* density, double* importance) {
    ErrorCode err;
    DagMC* dag = DAG;

    // the defaults from DagMC's old get_volume_metadata: mat = 0, rho = 0, imp = 1
    int mat_id = 0;
    double rho = 0, imp = 1;

    std::string str;
    err = dag->prop_value(vol, "mat", str);
    if(err == moab::MB_SUCCESS) {
        mat_id = strtol(str.c_str(), NULL, 10);
    }
    else if(err != moab::MB_TAG_NOT_FOUND) {
        // TAG_NOT_FOUND should not be returned as an error; it just means
        // the default value of mat_id needs to be used.
        CHECKERR(err);
    }
    
    err = dag->prop_value(vol, "rho", str);
    if(err == moab::MB_SUCCESS) {
        rho = strtod(str.c_str(), NULL);
    }
    else if(err != moab::MB_TAG_NOT_FOUND) {
        CHECKERR(err);
    }

    err = dag->prop_value(vol, "imp", str);
    if(err == moab::MB_SUCCESS) {
        imp = strtod(str.c_str(), NULL);
    }
    else if(err != moab::MB_TAG_NOT_FOUND) {
        CHECKERR(err);
    }

    *material = mat_id;
    *density = rho;
    *importance = imp;

    return moab::MB_SUCCESS;
}

ErrorCode get_volume_boundary(EntityHandle vol, vec3 minPt, vec3 maxPt) {
    return DAG->getobb(vol, minPt, maxPt);
}

} // namespace pyne
