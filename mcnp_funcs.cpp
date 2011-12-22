#include "mcnp_funcs.h"

#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "DagMC.hpp"
using moab::DagMC;

#include <limits>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef CUBIT_LIBS_PRESENT
#include <fenv.h>
#endif

  // globals

#define DAG DagMC::instance()

#define DGFM_SEQ   0
#define DGFM_READ  1
#define DGFM_BCAST 2

#ifdef ENABLE_RAYSTAT_DUMPS

#include <fstream>
#include <numeric>

static std::ostream* raystat_dump = NULL;

#endif 


/* Static values used by dagmctrack_ */

static DagMC::RayHistory history;
static int last_nps = 0;
static double last_uvw[3] = {0,0,0};
static std::vector< DagMC::RayHistory > history_bank;
static std::vector< DagMC::RayHistory > pblcm_history_stack;
static bool visited_surface = false;

static bool use_dist_limit = false;
static double dist_limit; // needs to be thread-local


void dagmcinit_(char *cfile, int *clen,  // geom
                char *ftol,  int *ftlen, // faceting tolerance
                int *parallel_file_mode, // parallel read mode
                double* dagmc_version, int* moab_version, int* max_pbl )
{
 
  MBErrorCode rval;

#ifdef ENABLE_RAYSTAT_DUMPS
  // the file to which ray statistics dumps will be written
  raystat_dump = new std::ofstream("dagmc_raystat_dump.csv");
#endif 
  
  *dagmc_version = DAG->version();
  *moab_version = DAG->interface_revision();
  
    // terminate all filenames with null char
  cfile[*clen] = ftol[*ftlen] = '\0';

    // initialize this as -1 so that DAGMC internal defaults are preserved
    // user doesn't set this
  double arg_facet_tolerance = -1;
                                                                        
  if ( *ftlen > 0 ) arg_facet_tolerance = atof(ftol);
  
  // read geometry
  rval = DAG->load_file(cfile, arg_facet_tolerance );
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to read input file: " << cfile << std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef CUBIT_LIBS_PRESENT
  // The Cubit 10.2 libraries enable floating point exceptions.  
  // This is bad because MOAB may divide by zero and expect to continue executing.
  // See MOAB mailing list discussion on April 28 2010.
  // As a workaround, put a hold exceptions when Cubit is present.

  fenv_t old_fenv;
  if ( feholdexcept( &old_fenv ) ){
    std::cerr << "Warning: could not hold floating-point exceptions!" << std::endl;
  }
#endif

 
  // initialize geometry
  rval = DAG->init_OBBTree();
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  pblcm_history_stack.resize( *max_pbl+1 ); // fortran will index from 1

}

void dagmcwritefacets_(char *ffile, int *flen)  // facet file
{
    // terminate all filenames with null char
  ffile[*flen]  = '\0';

  MBErrorCode rval = DAG->write_mesh(ffile,*flen);
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to write mesh file: " << ffile <<  std::endl;
    exit(EXIT_FAILURE);
  }

  return;

}



void dagmcwritemcnp_(char *lfile, int *llen)  // file with cell/surface cards
                     
{
  MBErrorCode rval;

  lfile[*llen]  = '\0';

  // parse data from geometry
  rval = DAG->parse_metadata();
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  std::string lfname(lfile, *llen);
  // only overwrite mcnp log file if lfname==lcad, that is lfname.compare("lcad") == 0
  rval = DAG->write_mcnp(lfname, lfname.compare("lcad") == 0);
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to write mcnp file: " << lfile <<  std::endl;
    exit(EXIT_FAILURE);
  }

}

void dagmcangl_(int *jsu, double *xxx, double *yyy, double *zzz, double *ang)
{
  MBEntityHandle surf = DAG->entity_by_index( 2, *jsu );
  double xyz[3] = {*xxx, *yyy, *zzz};
  MBErrorCode rval = DAG->get_angle(surf, xyz, ang, &history );
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC: failed in calling get_angle" <<  std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef TRACE_DAGMC_CALLS
  std::cout << "angl: " << *xxx << ", " << *yyy << ", " << *zzz << " --> " 
            << ang[0] <<", " << ang[1] << ", " << ang[2] << std::endl;
  MBCartVect uvw(last_uvw);
  MBCartVect norm(ang);
  double aa = angle(uvw,norm) * (180.0/M_PI);
  std::cout << "    : " << aa << " deg to uvw" << (aa>90.0? " (!)":"")  << std::endl;
#endif
  
}

void dagmcchkcel_by_angle_( double *uuu, double *vvv, double *www, 
                            double *xxx, double *yyy, double *zzz,
                            int *jsu, int *i1, int *j)
{


#ifdef TRACE_DAGMC_CALLS
  std::cout<< " " << std::endl;
  std::cout<< "chkcel_by_angle: vol=" << DAG->id_by_index(3,*i1) << " surf=" << DAG->id_by_index(2,*jsu)
           << " xyz=" << *xxx  << " " << *yyy << " " << *zzz << std::endl;
  std::cout<< "               : uvw = " << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

  double xyz[3] = {*xxx, *yyy, *zzz};
  double uvw[3] = {*uuu, *vvv, *www};

  MBEntityHandle surf = DAG->entity_by_index( 2, *jsu );
  MBEntityHandle vol  = DAG->entity_by_index( 3, *i1 );

  int result;
  MBErrorCode rval = DAG->test_volume_boundary( vol, surf, xyz, uvw, result, &history );
  if( MB_SUCCESS != rval ){
    std::cerr << "DAGMC: failed calling test_volume_boundary" << std::endl;
    exit(EXIT_FAILURE);
  }

  switch (result)
      {
      case 1: 
        *j = 0; // inside==  1 -> inside volume -> j=0
        break;
      case 0:
        *j = 1; // inside== 0  -> outside volume -> j=1
        break;
      default:
        std::cerr << "Impossible result in dagmcchkcel_by_angle" << std::endl;
        exit(EXIT_FAILURE);
      }
 
#ifdef TRACE_DAGMC_CALLS
  std::cout<< "chkcel_by_angle: j=" << *j << std::endl;
#endif

}

void dagmcchkcel_(double *uuu,double *vvv,double *www,double *xxx,
                  double *yyy,double *zzz, int *i1, int *j)
{


#ifdef TRACE_DAGMC_CALLS
  std::cout<< " " << std::endl;
  std::cout<< "chkcel: vol=" << DAG->id_by_index(3,*i1) << " xyz=" << *xxx 
           << " " << *yyy << " " << *zzz << std::endl;
  std::cout<< "      : uvw = " << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

  int inside;
  MBEntityHandle vol = DAG->entity_by_index( 3, *i1 );
  double xyz[3] = {*xxx, *yyy, *zzz};
  double uvw[3] = {*uuu, *vvv, *www};
  MBErrorCode rval = DAG->point_in_volume( vol, xyz, inside, uvw );

  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC: failed in point_in_volume" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  if (MB_SUCCESS != rval) *j = -2;
  else
    switch (inside)
      {
      case 1: 
        *j = 0; // inside==  1 -> inside volume -> j=0
        break;
      case 0:
        *j = 1; // inside== 0  -> outside volume -> j=1
        break;
      case -1:
        *j = 1; // inside== -1 -> on boundary -> j=1 (assume leaving volume)
        break;
      default:
        std::cerr << "Impossible result in dagmcchkcel" << std::endl;
        exit(EXIT_FAILURE);
      }
  
#ifdef TRACE_DAGMC_CALLS
  std::cout<< "chkcel: j=" << *j << std::endl;
#endif

}


void dagmcdbmin_( int *ih, double *xxx, double *yyy, double *zzz, double *huge, double* dbmin)
{
  double point[3] = {*xxx, *yyy, *zzz};

  // get handle for this volume (*ih)
  MBEntityHandle vol  = DAG->entity_by_index( 3, *ih );

  // get distance to closest surface
  MBErrorCode rval = DAG->closest_to_location(vol,point,*dbmin);

  // if failed, return 'huge'
  if (MB_SUCCESS != rval) {
    *dbmin = *huge;
    std::cerr << "DAGMC: error in closest_to_location, returning huge value from dbmin_" <<  std::endl;
  }

#ifdef TRACE_DAGMC_CALLS
  std::cout << "dbmin " << DAG->id_by_index( 3, *ih ) << " dist = " << *dbmin << std::endl;
#endif

}

void dagmcnewcel_( int *jsu, int *icl, int *iap )
{

  MBEntityHandle surf = DAG->entity_by_index( 2, *jsu );
  MBEntityHandle vol  = DAG->entity_by_index( 3, *icl );
  MBEntityHandle newvol = 0;

  MBErrorCode rval = DAG->next_vol( surf, vol, newvol );
  if( MB_SUCCESS != rval ){
    *iap = -1;
    std::cerr << "DAGMC: error calling next_vol, newcel_ returning -1" << std::endl;
  }
  
  *iap = DAG->index_by_handle( newvol );

  visited_surface = true;
  
#ifdef TRACE_DAGMC_CALLS
  std::cout<< "newcel: prev_vol=" << DAG->id_by_index(3,*icl) << " surf= " 
           << DAG->id_by_index(2,*jsu) << " next_vol= " << DAG->id_by_index(3,*iap) <<std::endl;

#endif
}

void dagmc_surf_reflection_( double *uuu, double *vvv, double *www, int* verify_dir_change )
{


#ifdef TRACE_DAGMC_CALLS
  // compute and report the angle between old and new
  MBCartVect oldv(last_uvw);
  MBCartVect newv( *uuu, *vvv, *www );
  
  std::cout << "surf_reflection: " << angle(oldv,newv)*(180.0/M_PI) << std::endl;;
#endif

  // a surface was visited
  visited_surface = true;

  bool update = true;
  if( *verify_dir_change ){
    if( last_uvw[0] == *uuu && last_uvw[1] == *vvv && last_uvw[2] == *www  )
      update = false;
  }

  if( update ){
    last_uvw[0] = *uuu;
    last_uvw[1] = *vvv;
    last_uvw[2] = *www;
    history.reset_to_last_intersection();  
  }

#ifdef TRACE_DAGMC_CALLS
  else{
    // mark it in the log if nothing happened
    std::cout << "(noop)";
  }

  std::cout << std::endl;
#endif 

}

void dagmc_particle_terminate_( )
{
  history.reset();

#ifdef TRACE_DAGMC_CALLS
  std::cout << "particle_terminate:" << std::endl;
#endif
}

// *ih              - volue index
// *uuu, *vvv, *www - ray direction
// *xxx, *yyy, *zzz - ray point
// *huge            - passed to ray_fire as 'huge'
// *dls             - output from ray_fire as 'dist_traveled'
// *jap             - intersected surface index, or zero if none
// *jsu             - previous surface index
void dagmctrack_(int *ih, double *uuu,double *vvv,double *www,double *xxx,
                 double *yyy,double *zzz,double *huge,double *dls,int *jap,int *jsu,
                 int *nps )
{
    // Get data from IDs
  MBEntityHandle vol = DAG->entity_by_index( 3, *ih );
  MBEntityHandle prev = DAG->entity_by_index( 2, *jsu );
  MBEntityHandle next_surf = 0;
  double next_surf_dist;

#ifdef ENABLE_RAYSTAT_DUMPS
  moab::OrientedBoxTreeTool::TrvStats trv;
#endif 

  double point[3] = {*xxx,*yyy,*zzz};
  double dir[3]   = {*uuu,*vvv,*www};  

  /* detect streaming or reflecting situations */
  if( last_nps != *nps || prev == 0 ){
    // not streaming or reflecting: reset history
    history.reset(); 
#ifdef TRACE_DAGMC_CALLS
    std::cout << "track: new history" << std::endl;
#endif

  }
  else if( last_uvw[0] == *uuu && last_uvw[1] == *vvv && last_uvw[2] == *www ){
    // streaming -- use history without change 
    // unless a surface was not visited
    if( !visited_surface ){ 
      history.rollback_last_intersection();
#ifdef TRACE_DAGMC_CALLS
      std::cout << "     : (rbl)" << std::endl;
#endif
    }
#ifdef TRACE_DAGMC_CALLS
    std::cout << "track: streaming " << history.size() << std::endl;
#endif
  }
  else{
    // not streaming or reflecting
    history.reset();

#ifdef TRACE_DAGMC_CALLS
    std::cout << "track: reset" << std::endl;
#endif

  }

  MBErrorCode result = DAG->ray_fire(vol, point, dir, 
                                     next_surf, next_surf_dist, &history, 
                                     (use_dist_limit ? dist_limit : 0 )
#ifdef ENABLE_RAYSTAT_DUMPS
                                     , raystat_dump ? &trv : NULL 
#endif
                                     );

  
  if(MB_SUCCESS != result){
    std::cerr << "DAGMC: failed in ray_fire" << std::endl;
    exit( EXIT_FAILURE );
  }

  
  for( int i = 0; i < 3; ++i ){ last_uvw[i] = dir[i]; } 
  last_nps = *nps;

  // Return results: if next_surf exists, then next_surf_dist will be nearer than dist_limit (if any)
  if( next_surf != 0 ){
    *jap = DAG->index_by_handle( next_surf ); 
    *dls = next_surf_dist; 
  }
  else{
    // no next surface 
    *jap = 0;
    if( use_dist_limit ){
      // Dist limit on: return a number bigger than dist_limit
      *dls = dist_limit * 2.0;
    }
    else{
      // Dist limit off: return huge value, triggering lost particle
      *dls = *huge;
    }
  }

  visited_surface = false;
  
#ifdef ENABLE_RAYSTAT_DUMPS
  if( raystat_dump ){

    *raystat_dump << *ih << ",";
    *raystat_dump << trv.ray_tri_tests() << ",";
    *raystat_dump << std::accumulate( trv.nodes_visited().begin(), trv.nodes_visited().end(), 0 ) << ",";
    *raystat_dump << std::accumulate( trv.leaves_visited().begin(), trv.leaves_visited().end(), 0 ) << std::endl;

  }
#endif 

#ifdef TRACE_DAGMC_CALLS

  std::cout<< "track: vol=" << DAG->id_by_index(3,*ih) << " prev_surf=" << DAG->id_by_index(2,*jsu) 
           << " next_surf=" << DAG->id_by_index(2,*jap) << " nps=" << *nps <<std::endl;
  std::cout<< "     : xyz=" << *xxx << " " << *yyy << " "<< *zzz << " dist = " << *dls << std::flush;
  if( use_dist_limit && *jap == 0 ) std::cout << " > distlimit" << std::flush;
  std::cout << std::endl;
  std::cout<< "     : uvw=" << *uuu << " " << *vvv << " "<< *www << std::endl;
#endif

}

void dagmc_bank_push_( int* nbnk )
{
  if( ((unsigned)*nbnk) != history_bank.size() ){
    std::cerr << "bank push size mismatch: F" << *nbnk << " C" << history_bank.size() << std::endl;
  }
  history_bank.push_back( history );

#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_push (" << *nbnk+1 << ")" << std::endl;
#endif
}

void dagmc_bank_usetop_( ) 
{

#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_usetop" << std::endl;
#endif

  if( history_bank.size() ){
    history = history_bank.back();
  }
  else{
    std::cerr << "dagmc_bank_usetop_() called without bank history!" << std::endl;
  }
}

void dagmc_bank_pop_( int* nbnk )
{

  if( ((unsigned)*nbnk) != history_bank.size() ){
    std::cerr << "bank pop size mismatch: F" << *nbnk << " C" << history_bank.size() << std::endl;
  }

  if( history_bank.size() ){
    history_bank.pop_back( ); 
  }

#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_pop (" << *nbnk-1 << ")" << std::endl;
#endif

}

void dagmc_bank_clear_( )
{
  history_bank.clear();
#ifdef TRACE_DAGMC_CALLS
  std::cout << "bank_clear" << std::endl;
#endif
}

void dagmc_savpar_( int* n )
{
#ifdef TRACE_DAGMC_CALLS
  std::cout << "savpar: " << *n << " ("<< history.size() << ")" << std::endl;
#endif
  pblcm_history_stack[*n] = history;
}

void dagmc_getpar_( int* n )
{
#ifdef TRACE_DAGMC_CALLS
  std::cout << "getpar: " << *n << " (" << pblcm_history_stack[*n].size() << ")" << std::endl;
#endif
  history = pblcm_history_stack[*n];
}


void dagmcvolume_(int* mxa, double* vols, int* mxj, double* aras)
{
  MBErrorCode rval;
  
    // get size of each volume
  int num_vols = DAG->num_entities(3);
  for (int i = 0; i < num_vols; ++i) {
    rval = DAG->measure_volume( DAG->entity_by_index(3, i+1), vols[i*2] );
    if( MB_SUCCESS != rval ){
      std::cerr << "DAGMC: could not measure volume " << i+1 << std::endl;
      exit( EXIT_FAILURE );
    }
  }
  
    // get size of each surface
  int num_surfs = DAG->num_entities(2);
  for (int i = 0; i < num_surfs; ++i) {
    rval = DAG->measure_area( DAG->entity_by_index(2, i+1), aras[i*2] );
    if( MB_SUCCESS != rval ){
      std::cerr << "DAGMC: could not measure surface " << i+1 << std::endl;
      exit( EXIT_FAILURE );
    }
  }

}

void dagmc_setdis_(double *d)
{
  dist_limit = *d;
#ifdef TRACE_DAGMC_CALLS
  std::cout << "setdis: " << *d << std::endl;
#endif
}

void dagmc_set_settings_(int* fort_use_dist_limit, int* use_cad, double* overlap_thickness, int* srccell_mode )
{

  if( *fort_use_dist_limit ){
    std::cout << "DAGMC distance limit optimization is ENABLED" << std::endl;
    use_dist_limit = true;
  }

  if( *srccell_mode ){
    std::cout << "DAGMC source cell optimization is ENABLED (warning: experimental!)" << std::endl;
  }

  DAG->set_use_CAD( *use_cad );

  DAG->set_overlap_thickness( *overlap_thickness );

}

void dagmc_init_settings_(int* fort_use_dist_limit, int* use_cad,    
                          double* overlap_thickness, double* facet_tol, int* srccell_mode )
{

  *fort_use_dist_limit = use_dist_limit ? 1 : 0;

  *use_cad = DAG->use_CAD() ? 1 : 0;

  *overlap_thickness = DAG->overlap_thickness();
  
  *facet_tol = DAG->faceting_tolerance();


  if( *srccell_mode ){
    std::cout << "DAGMC source cell optimization is ENABLED (warning: experimental!)" << std::endl;
  }

}

void dagmc_version_(double* dagmcVersion)
{
  *dagmcVersion = DAG->version();
}

