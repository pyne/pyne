#include "dagmc_bridge.h"

#include <DagMC.hpp>

#include <vector>

using moab::DagMC;
using moab::EntityHandle;

#define CHECKERR( err ) \
    if( (err) != moab::MB_SUCCESS ) return err;

float dag_version( void ){
    return DagMC::version();
}

unsigned dag_rev_version( void ){
    return DagMC::interface_revision();
}

int dag_ent_handle_size( void ){
    return sizeof( EntityHandle );
}

static std::vector<int> surfList;
static std::vector<int> volList;

const int* geom_id_list( int dimension, int* number_of_items ){
    switch(dimension){
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

EntityHandle handle_from_id( int dimension, int id ){
    return DagMC::instance()->entity_by_id( dimension, id );
}

int id_from_handle( EntityHandle eh ){
    return DagMC::instance()->get_entity_id( eh );
}

ErrorCode dag_load( const char* filename ){
    ErrorCode err;

    DagMC* dag = DagMC::instance();

    err = dag->load_file( filename );
    CHECKERR( err );
    err = dag->init_OBBTree();
    CHECKERR( err );
    err = dag->parse_metadata();
    CHECKERR( err );

    int num_surfs = dag->num_entities( 2 );
    surfList.reserve( num_surfs );
    for( int i = 1; i <= num_surfs; ++i ){
        surfList.push_back( dag->id_by_index( 2, i ) );
    }

    int num_vols = dag->num_entities( 3 );
    volList.reserve( num_vols );
    for( int i = 1; i <= num_vols; ++i ){
        volList.push_back( dag->id_by_index( 3, i ) );
    }

    return err;
}


void* dag_alloc_ray_history( void ){
    return new DagMC::RayHistory();
}

void dag_dealloc_ray_history( void* r ){
    delete ( static_cast<DagMC::RayHistory*>(r) );
}

ErrorCode dag_ray_fire( EntityHandle vol, vec3 ray_start, vec3 ray_dir, EntityHandle* next_surf_ent, double* next_surf_dist,
                        void* history, double distance_limit ){
    ErrorCode err;

    DagMC*  dag = DagMC::instance();

    err = dag->ray_fire( vol, ray_start, ray_dir, *next_surf_ent, *next_surf_dist, 
                         static_cast<DagMC::RayHistory*>(history), distance_limit );
    CHECKERR( err );

    return err;
}

ErrorCode dag_pt_in_vol( EntityHandle vol, vec3 pt, int* result, vec3 dir, const void* history){
    
    ErrorCode err;

    DagMC* dag = DagMC::instance();
    
    err = dag->point_in_volume( vol, pt, *result, dir, static_cast<const DagMC::RayHistory*>(history) );

    return err;
}

ErrorCode dag_next_vol( EntityHandle surface, EntityHandle volume, EntityHandle* next_vol ){

    ErrorCode err;
    DagMC* dag = DagMC::instance();

    err = dag->next_vol( surface, volume, *next_vol );

    return err;
}

int vol_is_graveyard( EntityHandle vol ){
    return DagMC::instance()->is_graveyard( vol );
}

/* int surf_is_spec_refl( EntityHandle surf ); */
/* int surf_is_white_refl( EntityHandle surf ); */

int vol_is_implicit_complement( EntityHandle vol ){
    return DagMC::instance()->is_implicit_complement( vol );
}

ErrorCode get_volume_metadata( EntityHandle vol, int* material, double* density, double* importance ){
    DagmcVolData vd;
    ErrorCode err;
    DagMC* dag = DagMC::instance();

    err = dag->get_volume_metadata( vol, vd );

    CHECKERR( err );

    *material = vd.mat_id;
    *density = vd.density;
    *importance = vd.importance;

    return err;
}
