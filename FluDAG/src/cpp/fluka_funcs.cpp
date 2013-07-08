//----------------------------------*-C++, Fortran-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/fluka_funcs.cpp
 * \author Julie Zachman 
 * \date   Mon Mar 22 2013 
 * \brief  Functions called by fluka
 * \note   After mcnp_funcs
 */
//---------------------------------------------------------------------------//
// $Id: 
//---------------------------------------------------------------------------//

#include "fluka_funcs.h"
#include "fludag_utils.h"

#include "DagWrappers.hh"
#include "dagmc_utils.hpp"

#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "DagMC.hpp"
#include "moab/Types.hpp"
using moab::DagMC;

#include <limits>
#include <ios>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <set>
#include <cstring>

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
#define DEBUG 1
/* These 37 strings are predefined FLUKA materials. Any ASSIGNMAt of unique 
 * materials not on this list requires a MATERIAL card. */
std::string flukaMatStrings[] = {"BLCKHOLE", "VACUUM", "HYDROGEN",
"HELIUM", "BERYLLIU", "CARBON", "NITROGEN", "OXYGEN", "MAGNESIU",      
"ALUMINUM", "IRON", "COPPER", "SILVER", "SILICON", "GOLD", "MERCURY",  
"LEAD", "TANTALUM", "SODIUM", "ARGON", "CALCIUM", "TIN", "TUNGSTEN",   
"TITANIUM", "NICKEL", "WATER", "POLYSTYR", "PLASCINT", "PMMA",         
"BONECOMP", "BONECORT", "MUSCLESK", "MUSCLEST", "ADTISSUE", "KAPTON",  
"POLYETHY", "AIR"};

int NUM_FLUKA_MATS = 37;

/* Create a set out of the hardcoded string array. */
std::set<std::string> FLUKA_mat_set(flukaMatStrings, flukaMatStrings+NUM_FLUKA_MATS); 

/* Maximum character-length of a cubit-named material property */
int MAX_MATERIAL_NAME_SIZE = 32;

bool debug = false; //true ;

/* Static values used by dagmctrack_ */

static DagMC::RayHistory history;
static int last_nps = 0;
static double last_uvw[3] = {0,0,0};
static std::vector< DagMC::RayHistory > history_bank;
static std::vector< DagMC::RayHistory > pblcm_history_stack;
static bool visited_surface = false;

static bool use_dist_limit = false;
static double dist_limit; // needs to be thread-local

MBEntityHandle next_surf; // the next suface the ray will hit
MBEntityHandle prev_surf; // the last value of next surface
MBEntityHandle PrevRegion; // the integer region that the particle was in previously

std::string ExePath() 
{
    int MAX_PATH = 256;
    char buffer[MAX_PATH];
    getcwd(  buffer, MAX_PATH );
    // std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
    // return std::string( buffer ).substr( 0, pos);
    return std::string( buffer );
}
/**	
  dagmcinit_ is meant to be called from a fortran caller.  Strings have to be 
  accompanied by their length, and will need to be null-appended.
  Precondition:  cfile exists and is readable
*/
void dagmcinit_(char *cfile, int *clen,  // geom
                char *ftol,  int *ftlen, // faceting tolerance
                int *parallel_file_mode, // parallel read mode
                double* dagmc_version, int* moab_version, int* max_pbl )
{
        // Presumably this serves as output to a calling fortran program
        *dagmc_version= DAG->version();
        *moab_version = DAG->interface_revision();
        // terminate all filenames with null char
        cfile[*clen] = ftol[*ftlen] = '\0';
        // Call as if running with fluka (final param=true)
	cpp_dagmcinit(cfile, *parallel_file_mode, *max_pbl);
}


/** 
   cpp_dagmcinit is called directly from c++ or from a fortran-called wrapper.
  Precondition:  myfile exists and is readable
*/
void cpp_dagmcinit(std::string infile,         // geom
                int parallel_file_mode, // parallel read mode
                int max_pbl)
{
 
  MBErrorCode rval;

  // initialize this as -1 so that DAGMC internal defaults are preserved
  // user doesn't set this
  double arg_facet_tolerance = -1;
                                                                        
  // jcz: leave arg_facet_tolerance as defined on previous line
  // if ( *ftlen > 0 ) arg_facet_tolerance = atof(ftol);

  // read geometry
  std::cerr << "Loading " << infile << std::endl;
  rval = DAG->load_file(infile.c_str(), arg_facet_tolerance );
  if (MB_SUCCESS != rval) 
    {
      std::cerr << "DAGMC failed to read input file: " << infile << std::endl;
      exit(EXIT_FAILURE);
    }
 
  // initialize geometry
  rval = DAG->init_OBBTree();
  if (MB_SUCCESS != rval) 
    {
      std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
      exit(EXIT_FAILURE);
    }
}

/**************************************************************************************************/
/******                                FLUKA stubs                                         ********/
/**************************************************************************************************/
/// From Flugg Wrappers
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// jomiwr(..)
//---------------------------------------------------------------------------//
/// Initialization routine, was in WrapInit.c
void jomiwr(int & nge, const int& lin, const int& lou, int& flukaReg)
{
  if(debug)
    {
      std::cout << "================== JOMIWR =================" << std::endl;
    }
  // return FLUGG code to fluka
  // nge = 3;

  //Original comment:  returns number of volumes + 1
  unsigned int numVol = DAG->num_entities(3);
  flukaReg = numVol;

  if(debug)
    {
      std::cout << "Number of volumes: " << flukaReg << std::endl;
      std::cout << "================== Out of JOMIWR =================" << std::endl;
    }

  return;
}

//---------------------------------------------------------------------------//
// g_step(..)
//---------------------------------------------------------------------------//
//  returns approved step of particle and all variables 
//
void g_step(double& pSx, 
          double& pSy, 
          double& pSz, 
          double* pV,
          int& oldReg,         // pass through
          const int& oldLttc,  // ignore
          double& propStep,    // .
          int& nascFlag,       // .
          double& retStep,     // reset in this method
          int& newReg,         // return from callee
          double& saf,         // ignore
          int& newLttc,        // .
          int& LttcFlag,       // . 
          double* sLt,         // .
          int* jrLt)           // .
{
  double safety; // safety parameter

  if(debug)
    {
      std::cout<<"============= G_STEP	 =============="<<std::endl;    
      std::cout << "Position " << pSx << " " << pSy << " " << pSz << std::endl;
      std::cout << "Direction vector " << pV[0] << " " << pV[1] << " " << pV[2] << std::endl;
      std::cout << "Oldreg = " << oldReg << std::endl;
      std::cout << "PropStep = " << propStep << std::endl;
    }
  
  double point[3] = {pSx,pSy,pSz};
  double dir[3]   = {pV[0],pV[1],pV[2]};  

  g_fire(oldReg, point, dir, propStep, retStep, newReg); // fire a ray 

  /*
  std::cout << oldReg << " " << newReg << " " << std::endl;
  std::cout.precision(25);
  std::cout << std::scientific ;
  std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;
  std::cout << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
  std::cout << retStep << " " << propStep << std::endl;
  */

  if(debug)
    {
      std::cout << "saf = " << saf << std::endl;
      std::cout << std::setw(20) << std::scientific;
      std::cout << "newReg = " << newReg << " retStep = " << retStep << std::endl;
    }

  return;
}

//---------------------------------------------------------------------------//
// void g_fire(int& oldRegion, double point[], double dir[], 
//              double &propStep, double& retStep,  int& newRegion)
//---------------------------------------------------------------------------//
// oldRegion - the region of the particle's current coordinates
// point     - the particle's coordinate location vector
// dir       - the direction vector of the particle's current path (ray)
// propStep  - ??
// retStep   - returned as the distance from the particle's current location, along its ray, to the next boundary
// newRegion - gotten from the value returned by DAG->next_vol
// newRegion is gotten from the volue returned by DAG->next_vol
void g_fire(int& oldRegion, double point[], double dir[], double &propStep, double& retStep,  int& newRegion)
{

  /*
  if (PrevRegion != oldRegion) // if the particle is not in the correct volume since could be a banked history
    {
      int dummy;
      int errFlg;
      f_look(point[0],point[1],point[2],dir,0,dummy,oldRegion,errFlg,dummy); // where is particle
    }
  */

  if(debug)
  {
      std::cout<<"============= g_fire =============="<<std::endl;    
      std::cout << "Point " << point[0] << " " << point[1] << " " << point[2] << std::endl;
      std::cout << "Direction vector " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
  }
  MBEntityHandle vol = DAG->entity_by_index(3,oldRegion);
  double next_surf_dist;
  MBEntityHandle newvol = 0;

  next_surf = prev_surf;

  // vol = check_reg(vol,point,dir); // check we are where we say we are
  oldRegion = DAG->index_by_handle(vol);
  // next_surf is a global
  MBErrorCode result = DAG->ray_fire(vol, point, dir, next_surf, next_surf_dist );
  if ( result != MB_SUCCESS )
    {
      std::cout << "DAG ray fire error" << std::endl;
      exit(0);
    }

  retStep = next_surf_dist; // the returned step length is the distance to next surf

  if ( next_surf == 0 ) // if next_surface is 0 then we are lost
    {
      std::cout << "!!! Lost Particle !!! " << std::endl;
      std::cout << "in region, " << oldRegion << " aka " << DAG->entity_by_index(3,oldRegion) << std::endl;  

      std::cout.precision(25);
      std::cout << std::scientific ; 
      std::cout << "position of particle " << point[0] << " " << point[1] << " " << point[2] << std::endl;
      std::cout << " traveling in direction " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
      std::cout << "!!! Lost Particle !!!" << std::endl;
      exit(0);
    }
  
  if ( propStep > retStep ) // will cross into next volume next step
    {
      MBErrorCode rval = DAG->next_vol(next_surf,vol,newvol);
      newRegion = DAG->index_by_handle(newvol);
      retStep = retStep ; // path limited by geometry
      if ( retStep < 1.0e-9 )
	retStep = 1.0e-9;
    }
  else
    {
      newRegion = oldRegion;
      retStep = propStep; //physics limits step
      next_surf = prev_surf;
    }

  PrevRegion = newRegion; // particle will be moving to PrevRegion upon next entry.

  if(debug)
  {
     std::cout << "Region on other side of surface is  = " << newRegion << \
                  ", Distance to next surf is " << retStep << std::endl;
  }

  prev_surf = next_surf;

  return;
}
///////			End g_step and g_fire
/////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------//
// normal
//---------------------------------------------------------------------------//
// Local wrapper for fortran-called, f_normal.  This function is supplied for testing
// purposes.  Its signature shows what parameters are being used in our wrapper 
// implementation.  
// Any FluDAG calls to f_normal should use this call instead.
// ASSUMES:  no ray history
// Notes
// - direction is not taken into account 
// - curRegion is not currently used.  It is expected to be implemented as a check
//   on what the sign of the normal should be.  It is used in a call to DAG->surface_sense
int  normal (double& posx, double& posy, double& posz, double *norml, int& curRegion)
{
   int flagErr; 
   int dummyReg;
   double dummyDirx, dummyDiry, dummyDirz;
   f_normal(posx, posy, posz, dummyDirx, dummyDiry, dummyDirz, norml, curRegion, dummyReg, flagErr);
   return flagErr;
}
//---------------------------------------------------------------------------//
// f_normal(..)
//---------------------------------------------------------------------------//
/// From Flugg Wrappers WrapNorml.cc
//  Note:  The normal is calculated at the point on the surface nearest the 
//         given point
// ASSUMES:  Point is on the boundary
// Parameters Set:
//     norml vector
//     flagErr = 0 if ok, !=0 otherwise
// Does NOT set any region, point or direction vector.
// Globals used:
//     next_surf, set by ray_fire 
void f_normal(double& pSx, double& pSy, double& pSz,
            double& pVx, double& pVy, double& pVz,
	    double* norml, const int& oldRegion, 
	    const int& newReg, int& flagErr)
{
  if(debug)
  {
      std::cout << "============ NRMLWR =============" << std::endl;
  }

  MBEntityHandle OldReg = DAG -> entity_by_index(3,oldRegion); // entity handle
  double xyz[3] = {pSx,pSy,pSz}; //position vector
  double uvw[3] = {pVx,pVy,pVz}; //particl directoin
  int result; // particle is entering or leaving

  MBErrorCode ErrorCode = DAG->test_volume_boundary( OldReg, next_surf,xyz,uvw, result);  // see if we are on boundary
  ErrorCode = DAG->get_angle(next_surf,xyz,norml); 
  // result = 1 entering, 0 leaving
  if ( result == 0 ) // vector should point towards OldReg
    {
      norml[0] = norml[0]*-1.0;
      norml[1] = norml[1]*-1.0;
      norml[2] = norml[2]*-1.0;
    }

  if(debug)
  {
      std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
  }
  return;
}
//---------------------------------------------------------------------------//
// getSense(..)
//---------------------------------------------------------------------------//
// Helper function
int getSense(int region)
{

  int sense;  // sense of next_surf with respect to oldRegion (volume)

  MBEntityHandle vol = DAG->entity_by_index(3, region);
 
  MBErrorCode ErrorCode = DAG->surface_sense(vol, next_surf, sense); 
  if(false)
  {
      std::cout << "Sense of next_surf with respect to the region is " << sense << std::endl;
  }
  return sense; 
} 
///////			End f_normal, normal, and getSense(..)
/////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------//
// look(..)
//---------------------------------------------------------------------------//
// Testable local wrapper for fortran-called, f_look
// This function signature shows what parameters are being used in our wrapper implementation
// oldRegion is looked at if we are no a boundary, but it is not set.
// ASSUMES:  position is not on a boundary
// RETURNS: nextRegion, the region the given point is in 
int look( double& posx, double& posy, double& posz, double* dir, int& oldRegion)
{
   int flagErr;
   int lattice_dummy;  // not used
   int nextRegion;     
   f_look(posx, posy, posz, dir, oldRegion, lattice_dummy, nextRegion, flagErr, lattice_dummy);
   return nextRegion;
}
//---------------------------------------------------------------------------//
// f_look(..)
//---------------------------------------------------------------------------//
// Wrapper for localisation of starting point of particle.
//
// Question:  Should pV, the direction vector, be used?  The Flugg wrapper
//            code passes it to a geant call 
//           "ptrNavig->LocateGlobalPointAndUpdateTouchable()"
//////////////////////////////////////////////////////////////////
// This function answers the question What volume is the point in?  
// oldReg - not looked at UNLESS the volume is on the boundary, then newReg=oldReg
// nextRegion - set to the volume index the point is in.
// ToDo:  Is there an error condition for the flagErr that is guaranteed not to be equal to the next region?
//        Find a way to make use of the error return from point_in_volume
void f_look(double& pSx, double& pSy, double& pSz,
          double* pV, const int& oldReg, const int& oldLttc,
          int& nextRegion, int& flagErr, int& newLttc)
{
  if(debug)
  {
      std::cout << "======= LKWR =======" << std::endl;
      std::cout << "position is " << pSx << " " << pSy << " " << pSz << std::endl; 
  }

  double xyz[] = {pSx, pSy, pSz};       // location of the particle (xyz)
  const double dir[] = {pV[0],pV[1],pV[2]};
  // Initialize to outside boundary.  This value can be 0 or +/-1 for ouside, inside, or on boundary.
  // ToDo:  Should this be initialized at all?  Or should it be initialized to an invalide number?
  int is_inside = 0;                    
  int num_vols = DAG->num_entities(3);  // number of volumes

  for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
    {
      MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
      // No ray history 
      MBErrorCode code = DAG->point_in_volume(volume, xyz, is_inside,dir);

      // check for non error
      if(MB_SUCCESS != code) 
	{
	  std::cout << "Error return from point_in_volume!" << std::endl;
	  flagErr = -33;
	  return;
	}

      if ( is_inside == 1 ) // we are inside the cell tested
	{
	  nextRegion = i;
          //BIZARRELY - WHEN WE ARE INSIDE A VOLUME, BOTH, nextRegion has to equal flagErr
	  flagErr = nextRegion;
	  return;	  
	}
      else if ( is_inside == -1 )
	{
	  std::cout << "We cannot be here" << std::endl;
	  exit(0);
	}
    }  // end loop over all volumes

  // ToDo:  remove?  jcz
  //special_check(xyz,dir,nextRegion);
  // if we return update xyz
  //pSx=xyz[0];
  //pSy=xyz[1];
  //pSz=xyz[2];
  // if we are here do slow check
  slow_check(xyz,dir,nextRegion);
  flagErr = nextRegion; // return nextRegion
  return;
}

//---------------------------------------------------------------------------//
// slow_check(..)
//---------------------------------------------------------------------------//
// Helper function
void slow_check(double pos[3], const double dir[3], int &oldReg)
{
  std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
  std::cout << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
  int num_vols = DAG->num_entities(3);  // number of volumes
  int is_inside = 0;
  for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
    {
      MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
      MBErrorCode code = DAG->point_in_volume(volume, pos, is_inside,dir); 
      if ( code != MB_SUCCESS)
	{
	 std::cout << "Failure from point in volume" << std::endl;
	 exit(0);
	}

      if ( is_inside == 1) // if in volume
	{
	  oldReg = DAG->index_by_handle(volume); //set oldReg
	  std::cout << pos[0] << " " << pos[1] << " " << pos[2] << " " << oldReg << std::endl;
	  return;
	}
    }

  std::cout << "FAILED SLOW CHECK" << std::endl;
  exit(0);
}

//---------------------------------------------------------------------------//
// check_reg(..)
//---------------------------------------------------------------------------//
// NOT CALLED - Helper function
// check we are where we say we are
MBEntityHandle check_reg(MBEntityHandle volume, double point[3], double dir[3]) 
{
  int is_inside;
  MBErrorCode code = DAG->point_in_volume(volume, point, is_inside,dir); 
  if (is_inside == 1 )
    {
      // we are where we say we are
      return volume;
    }
  else
    {
      int num_vols = DAG->num_entities(3);  // number of volumes
      for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
	{
	  MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
	  MBErrorCode code = DAG->point_in_volume(volume, point, is_inside,dir); 
	  if ( is_inside == 1) // if in volume
	    {
	      return volume;
	    }
	}
      std::cout.precision(25);
      std::cout << std::scientific ; 
      std::cout << "position of particle " << point[0] << " " << point[1] << " " << point[2] << std::endl;
      std::cout << " traveling in direction " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
      std::cout << "particle not nowhere" << std::endl;
      exit(0);
    }
}

// ToDo:  remove? jcz
//---------------------------------------------------------------------------//
// special_check(..)
//---------------------------------------------------------------------------//
// NOT CALLED - Helper function
/*
void special_check(double pos[3],const double dir[3], int& oldReg)
{
  int num_vols = DAG->num_entities(3);  // number of volumes
  int counter = 0; //
  int is_inside = 0;
  do 
    {
      // bump particle position along dir
      pos[0]=pos[0]+(dir[0]*1.0e-9);
      pos[1]=pos[1]+(dir[1]*1.0e-9);
      pos[2]=pos[2]+(dir[2]*1.0e-9);

      for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
	{
	  MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
	  MBErrorCode code = DAG->point_in_volume(volume, pos, is_inside,dir); 
	  if ( is_inside == 1) // if in volume
	    {
	      std::cout << "had to bump " << counter << " times" << std::endl;
	      oldReg = DAG->index_by_handle(volume); //set oldReg
	      return;
	    }
	}
      counter++;
    }
  while ( is_inside != 0 );
}
*/
void f_lookdb(double& pSx, double& pSy, double& pSz,
	    double* pV, const int& oldReg, const int& oldLttc,
	    int& newReg, int& flagErr, int& newLttc)
{
  if(debug)
    {
      std::cout<<"============= F_LooKDB =============="<< std::endl;
    }
  //return region number and dummy variables
  newReg=0;   
  newLttc=0;
  flagErr=-1; 

  return;
}


/*
 * f_g1rt
 */
void f_g1rt(void)
{
  if(debug)
    {
      std::cout<<"============ F_G1RT ============="<<std::endl;
    }
    return;
}

// Set DNEAR option if needed
int f_idnr(const int & nreg, const int & mlat) 

{
	

// returns 0 if user doesn't want Fluka to use DNEAR to compute the 
// step (the same effect is obtained with the GLOBAL (WHAT(3)=-1)
// card in fluka input), returns 1 if user wants Fluka always to use DNEAR.

	return 0;
}
/**************************************************************************************************/
/******                                End of FLUKA stubs                                  ********/
/**************************************************************************************************/

/*
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
*/
/**
 * Helper function for parsing DagMC properties that are integers.
 * Returns true on success, false if property does not exist on the volume,
 * in which case the result is unmodified.
 * If DagMC throws an error, calls exit().
 */
static bool get_int_prop( MBEntityHandle vol, int cell_id, const std::string& property, int& result ){

  MBErrorCode rval;
  if( DAG->has_prop( vol, property ) ){
    std::string propval;
    rval = DAG->prop_value( vol, property, propval );
    if( MB_SUCCESS != rval ){ 
      std::cerr << "DagMC failed to get expected property " << property << " on cell " << cell_id << std::endl;
      std::cerr << "Error code: " << rval << std::endl;
      exit( EXIT_FAILURE );
    }
    const char* valst = propval.c_str();
    char* valend;
    result = strtol( valst, &valend, 10 );
    if( valend[0] != '\0' ){
      // strtol did not consume whole string
      std::cerr << "DagMC: trouble parsing '" << property <<"' value (" << propval << ") for cell " << cell_id << std::endl;
      std::cerr << "       the parsed value is " << result << ", using that." << std::endl;
    }
    return true;
  }
  else return false;

}

/**
 * Helper function for parsing DagMC properties that are doubles.
 * Returns true on success, false if property does not exist on the volume,
 * in which case the result is unmodified.
 * If DagMC throws an error, calls exit().
 */
static bool get_real_prop( MBEntityHandle vol, int cell_id, const std::string& property, double& result ){

  MBErrorCode rval;
  if( DAG->has_prop( vol, property ) ){
    std::string propval;
    rval = DAG->prop_value( vol, property, propval );
    if( MB_SUCCESS != rval ){ 
      std::cerr << "DagMC failed to get expected property " << property << " on cell " << cell_id << std::endl;
      std::cerr << "Error code: " << rval << std::endl;
      exit( EXIT_FAILURE );
    }
    const char* valst = propval.c_str();
    char* valend;
    result = strtod( valst, &valend );
    if( valend[0] != '\0' ){
      // strtod did not consume whole string
      std::cerr << "DagMC: trouble parsing '" << property <<"' value (" << propval << ") for cell " << cell_id << std::endl;
      std::cerr << "       the parsed value is " << result << ", using that." << std::endl;
    }
    return true;
  }
  else return false;

}

// take a string like "surf.flux.n", a key like "surf.flux", and a number like 2,
// If the first part of the string matches the key, remove the key from the string (leaving, e.g. ".n")
// and return the number.
static int tallytype( std::string& str, const char* key, int ret )
{
  if( str.find(key) == 0 ){
    str.erase( 0, strlen(key) );
    return ret;
  }
  return 0;
}

//---------------------------------------------------------------------------//
// get_tallyspec(..)
//---------------------------------------------------------------------------//
/// Needed by fludagwrite
// given a tally specifier like "1.surf.flux.n", return a printable card for the specifier
// and set 'dim' to 2 or 3 depending on whether its a surf or volume tally
static char* get_tallyspec( std::string spec, int& dim ){

  if( spec.length() < 2 ) return NULL;
  const char* str = spec.c_str();
  char* p;

  int ID = strtol( str, &p, 10 );
  if( p == str ) return NULL; // did not find a number at the beginning of the string
  if( *p != '.' ) return NULL; // did not find required separator
  str = p + 1;

  if( strlen(str) < 1 ) return NULL;

  std::string tmod;
  if( str[0] == 'q' ){
    tmod = "+"; 
    str++;
  }
  else if( str[0] == 'e' ){
    tmod = "*";
    str++;
  }
  
  std::string remainder(str);
  int type = 0;
  type = tallytype( remainder, "surf.current", 1 );
  if(!type) type = tallytype( remainder, "surf.flux", 2 );
  if(!type) type = tallytype( remainder, "cell.flux", 4 );
  if(!type) type = tallytype( remainder, "cell.heating", 6 );
  if(!type) type = tallytype( remainder, "cell.fission", 7 );
  if(!type) type = tallytype( remainder, "pulse.height", 8 );
  if( type == 0 ) return NULL;

  std::string particle = "n";
  if( remainder.length() >= 2 ){
    if(remainder[0] != '.') return NULL;
    particle = remainder.substr(1);
  }
  
  char* ret = new char[80];
  sprintf( ret, "%sf%d:%s", tmod.c_str(), (10*ID+type), particle.c_str() );

  dim = 3;
  if( type == 1 || type == 2 ) dim = 2;
  return ret;

}

//---------------------------------------------------------------------------//
// fludagwrite_assignma
//---------------------------------------------------------------------------//
/// Called from mainFludag when only one argument is given to the program.
//  This function writes out a simple numerical material assignment to the named argument file
//  Example usage:  mainFludag dagmc.html
//  Outputs
//           mat.inp  contains MATERIAL and ASSIGNMAt records for the input geometry.
//                    The MATERIAL is gotten by parsing the Cubit volume name on underscores.  
//                    The string after "M_" is considered to be t he material for that volume.
//                    There are no MATERIAL cards for the materials in the FLUKA_mat_set list
//                    For the remaining materials, there is one MATERIAL card apiece (no dups)
//                    User-named (not predefined) materials are TRUNCATED to 8 chars.
//                    User-named material id's start at 25 and increment by 1 for each MATERIAL card
//           index-id.txt  Map of FluDAG volume index vs Cubit volume ids, for info only.
//  Note that a preprocessing step to this call sets up the the DAG object that contains 
//  all the geometry information contained in dagmc.html.  
//  the name of the (currently hardcoded) output file is "mat.inp"
//  The graveyard is assumed to be the last region.
void fludagwrite_assignma(std::string lfname)  // file with cell/surface cards
{
  int num_vols = DAG->num_entities(3);
  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::cout << "\tnum_vols is " << num_vols << std::endl;
  MBErrorCode ret;
  MBEntityHandle entity = 0;
  int id;

  std::vector< std::string > keywords;
  ret = DAG->detect_available_props( keywords );
  // parse data from geometry so that property can be found
  ret = DAG->parse_properties( keywords );

  if (MB_SUCCESS != ret) 
  {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }
  // Preprocessing loop:  make a string, "props",  for each vol entity
  // This loop could be removed if the user doesn't care to see terminal output
  std::cout << "Property list: " << std::endl;
  for (unsigned i=1; i<=num_vols; i++)
  {
     
      std::string props = mat_property_string(i, keywords);
      id = DAG->id_by_index(3, i);
      if (props.length()) 
      {
         std::cout << "Vol " << i << ", id=" << id << ": parsed props: " << props << std::endl; 
      }
      else
      {
         std::cout << "Vol " << i << ", id=" << id << " has no props: " <<  std::endl; 
      }
  }

  // Open an outputstring for mat.inp
  std::ostringstream ostr;
  // Open an outputstring for index-id table and put a header in it
  std::ostringstream idstr;
  idstr << std::setw(5) <<  "Index" ;
  idstr << std::setw(5) <<  "   Id" << std::endl;

  // Prepare a list to contain unique materials not in Flulka's list
  std::list<std::string> uniqueMatList;

  // Loop through 3d entities.  In model_complete.h5m there are 90 vols
  std::vector<std::string> vals;
  std::string material;
  char buffer[MAX_MATERIAL_NAME_SIZE];
  for (unsigned i = 1 ; i <= num_vols ; i++)
  {
      vals.clear();
      entity = DAG->entity_by_index(3, i);
      id = DAG->id_by_index(3, i);
      // Create the id-index string for this vol
      idstr << std::setw(5) << std::right << i;
      idstr << std::setw(5) << std::right << id << std::endl;
      // Create the mat.inp string for this vol
      if (DAG->has_prop(entity, "graveyard"))
      {
	 ostr << std::setw(10) << std::left  << "ASSIGNMAt";
	 ostr << std::setw(10) << std::right << "BLCKHOLE";
	 ostr << std::setw(10) << std::right << i << std::endl;
      }
      else if (DAG->has_prop(entity, "M"))
      {
         ret = DAG->prop_values(entity, "M", vals);
         if (vals.size() >= 1)
         {
            // Make a copy of string in vals[0]; full string needs to be compared to
            // FLUKA materials list; copy is for potential truncation
            std::strcpy(buffer, vals[0].c_str());
            material = std::string(buffer);
            
            if (vals[0].size() > 8)
            {
               material.resize(8);
            }
            if (FLUKA_mat_set.find(vals[0]) == FLUKA_mat_set.end())
            {
                // current material is not in the pre-existing FLUKA material list
                uniqueMatList.push_back(material); 
                std::cerr << "Adding material " << material << " to the MATERIAL card list" << std::endl;
            }
         }
         else
         {
            material = "moreThanOne";
         }
	 ostr << std::setw(10) << std::left  << "ASSIGNMAt";
	 ostr << std::setw(10) << std::right << material;
	 ostr << std::setw(10) << std::right << i << std::endl;
      }
  }
  // Finish the ostr with the implicit complement card
  std::string implicit_comp_comment = "* The next volume is the implicit complement";
  ostr << implicit_comp_comment << std::endl;
  ostr << std::setw(10) << std::left  << "ASSIGNMAt";
  ostr << std::setw(10) << std::right << "VACUUM";
  ostr << std::setw(10) << std::right << num_vols+1 << std::endl;

  // Process the uniqueMatList list so that it truly is unique
  uniqueMatList.sort();
  uniqueMatList.unique();
  // Print the final list
  if (debug)
  {
     std::list<std::string>::iterator it; 
     for (it=uniqueMatList.begin(); it!=uniqueMatList.end(); ++it)
     {
        std::cerr << *it << std::endl;
     }
  
     // Show the output string just created
     std::cout << ostr.str();
  }

  // Prepare an output file of the given name; put a header and the output string in it
  std::ofstream lcadfile( lfname.c_str());
  std::string header = "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...";
  if (uniqueMatList.size() != 0)
  {
     int matID = 25;
     lcadfile << header << std::endl;
     std::list<std::string>::iterator it; 
     for (it=uniqueMatList.begin(); it!=uniqueMatList.end(); ++it)
     {
        lcadfile << std::setw(10) << std::left << "MATERIAL";
        lcadfile << std::setw(10) << std::right << "";
        lcadfile << std::setw(10) << std::right << "";
        lcadfile << std::setw(10) << std::right << "";
        lcadfile << std::setw(10) << std::right << ++matID;
        lcadfile << std::setw(10) << std::right << "";
        lcadfile << std::setw(10) << std::right << "";
        lcadfile << std::setw(10) << std::left << *it << std::endl;
     }
  }
  lcadfile << header << std::endl;
  lcadfile << ostr.str();
  lcadfile.close();

  // Prepare an output file named "index_id.txt" for idstr
  std::string index_id_filename = "index_id.txt";
  std::ofstream index_id(index_id_filename.c_str());
  index_id << idstr.str();
  index_id.close(); 
  std::cerr << "Writing lcad file = " << lfname << std::endl; 
// Before opening file for writing, check for an existing file
/*
  if( lfname != "lcad" ){
    // Do not overwrite a lcad file if it already exists, except if it has the default name "lcad"
    if( access( lfname.c_str(), R_OK ) == 0 ){
      std::cout << "DagMC: reading from existing lcad file " << lfname << std::endl;
      return; 
    }
  }
*/

}

//---------------------------------------------------------------------------//
// mat_property_string
//---------------------------------------------------------------------------//
// For a given volume, find the values of all properties named "MAT".   
// Create a string with these properites
// Modified from make_property_string
// This function helps with debugging, but is not germane to writing cards.
std::string mat_property_string (int index, std::vector<std::string> &properties)
{
  ErrorCode ret;
  std::string propstring;
  EntityHandle entity = DAG->entity_by_index(3,index);
  int id = DAG->id_by_index(3, index);
  for (std::vector<std::string>::iterator p = properties.begin(); p != properties.end(); ++p)
  {
     if ( DAG->has_prop(entity, *p) )
     {
        std::vector<std::string> vals;
        ret = DAG->prop_values(entity, *p, vals);
        CHECKERR(*DAG, ret);
        propstring += *p;
        if (vals.size() == 1)
        {
 	   propstring += "=";
           propstring += vals[0];
        }
        else if (vals.size() > 1)
        {
 	   // this property has multiple values; list within brackets
           propstring += "=[";
	   for (std::vector<std::string>::iterator i = vals.begin(); i != vals.end(); ++i)
           {
	       propstring += *i;
               propstring += ",";
           }
           // replace the last trailing comma with a close bracket
           propstring[ propstring.length() -1 ] = ']';
        }
        propstring += ", ";
     }
  }
  if (propstring.length())
  {
     propstring.resize( propstring.length() - 2); // drop trailing comma
  }
  return propstring;
}

//---------------------------------------------------------------------------//
// make_property_string
//---------------------------------------------------------------------------//
// For a given volume, find all properties associated with it, and any and all 
//     values associated with each property
// Copied and modified from obb_analysis.cpp
static std::string make_property_string (EntityHandle eh, std::vector<std::string> &properties)
{
  ErrorCode ret;
  std::string propstring;
  for (std::vector<std::string>::iterator p = properties.begin(); p != properties.end(); ++p)
  {
     if ( DAG->has_prop(eh, *p) )
     {
        std::vector<std::string> vals;
        ret = DAG->prop_values(eh, *p, vals);
        CHECKERR(*DAG, ret);
        propstring += *p;
        if (vals.size() == 1)
        {
 	   propstring += "=";
           propstring += vals[0];
        }
        else if (vals.size() > 1)
        {
 	   // this property has multiple values; list within brackets
           propstring += "=[";
	   for (std::vector<std::string>::iterator i = vals.begin(); i != vals.end(); ++i)
           {
	       propstring += *i;
               propstring += ",";
           }
           // replace the last trailing comma with a close bracket
           propstring[ propstring.length() -1 ] = ']';
        }
        propstring += ", ";
     }
  }
  if (propstring.length())
  {
     propstring.resize( propstring.length() - 2); // drop trailing comma
  }
  return propstring;
}

//////////////////////////////////////////////////////////////////////////
/////////////
/////////////		region2name - modified from dagmcwrite
/////////////
//////////////////////////////////////////////////////////////////////////
void region2name(int volindex, char *vname )  // file with cell/surface cards
{
  MBErrorCode rval;

  std::vector< std::string > fluka_keywords;
  fluka_keywords.push_back( "mat" );
  fluka_keywords.push_back( "rho" );
  fluka_keywords.push_back( "comp" );
  fluka_keywords.push_back( "graveyard" );

  // parse data from geometry
  rval = DAG->parse_properties (fluka_keywords);
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

// std::ostringstream ostr;

  int cmat = 0;
  double crho;

  MBEntityHandle vol = DAG->entity_by_index( 3, volindex );
  int cellid = DAG->id_by_index( 3, volindex);

  bool graveyard = DAG->has_prop( vol, "graveyard" );

  std::ostringstream istr;
  if( graveyard )
  {
     istr << "BLCKHOLE";
     if( DAG->has_prop(vol, "comp") )
     {
       // material for the implicit complement has been specified.
       get_int_prop( vol, cellid, "mat", cmat );
       get_real_prop( vol, cellid, "rho", crho );
       std::cout << "Detected material and density specified for implicit complement: " << cmat <<", " << crho << std::endl;
     }
   }
   else if( DAG->is_implicit_complement(vol) )
   {
      istr << "mat_" << cmat;
      if( cmat != 0 ) istr << "_rho_" << crho;
   }
   else
   {
      int mat = 0;
      get_int_prop( vol, cellid, "mat", mat );

      if( mat == 0 )
      {
        istr << "0";
      }
      else
      {
        double rho = 1.0;
        get_real_prop( vol, cellid, "rho", rho );
        istr << "mat_" << mat << "_rho_" << rho;
      }
   }
   char *cstr = new char[istr.str().length()+1];
   std:strcpy(cstr,istr.str().c_str());
   vname = cstr;
}

//////////////////////////////////////////////////////////////////////////
/////////////
/////////////		fludagwrite_mat - modified from dagmcwrite
/////////////
//////////////////////////////////////////////////////////////////////////
void fludagwrite_mat(std::string fname )  // file with cell/surface cards
{
  MBErrorCode rval;

  std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  std::vector< std::string > fluka_keywords;

  fluka_keywords.push_back( "mat" );
  fluka_keywords.push_back( "comp" );
  fluka_keywords.push_back( "graveyard" );
  

  // parse data from geometry
  rval = DAG->parse_properties (fluka_keywords);
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  std::cerr << "Going to write file = " << fname << std::endl;
  // Before opening file for writing, check for an existing file
  if( fname != "lcad" ){
    // Do not overwrite a lcad file if it already exists, except if it has the default name "lcad"
    if( access( fname.c_str(), R_OK ) == 0 ){
      std::cout << "DagMC: reading from existing lcad file " << fname << std::endl;
      return; 
    }
  }

  int num_cells = DAG->num_entities( 3 );
  int cmat = 0;

  // Open an outputstring
  std::ostringstream ostr;

  // write the cell cards
  for( int i = 1; i <= num_cells; ++i ){

    MBEntityHandle vol = DAG->entity_by_index( 3, i );
    int cellid = DAG->id_by_index( 3, i );

    // for model_complete goes from 1-88, 93, 94
    ostr << std::setw(10) << std::left  << "ASSIGNMAT";

    bool graveyard = DAG->has_prop( vol, "graveyard" );

    if( graveyard )
    {
      // lcadfile << "BLCKHOLE";
      ostr << "BLCKHOLE";
      if( DAG->has_prop(vol, "comp") )
      {
        // material for the implicit complement has been specified.
        get_int_prop( vol, cellid, "mat", cmat );
      }
    }
    else if( DAG->is_implicit_complement(vol) )
    {
      // lcadfile << "mat_" << cmat;
      // lcadfile << " $ implicit complement";
      ostr << "mat_" << cmat;
      ostr << " $ implicit complement";
    }
    else
    {
      int mat = 0;
      get_int_prop( vol, cellid, "mat", mat );

      if( mat == 0 )
      {
        // lcadfile << "0";
        ostr << std::setw(10) << std::right << "0";
      }
      else
      {
        ostr << std::setw(10) << std::right << "mat_" << mat;
      }
    }

    ostr << std::setw(10) << std::right << cellid << std::endl;
  } // end iteration through cells

  // Show the output string just created
  std::cout << ostr.str();

  // Prepare an output file of the given name; put a header and the output string in it
  std::ofstream lcadfile( fname.c_str() );
  std::string header = "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7...";
  lcadfile << header << std::endl;
  lcadfile << ostr.str();

  std::cerr << "Writing lcad file = " << fname << std::endl;
}

//////////////////////////////////////////////////////////////////////////
/////////////
/////////////		fludagwrite modified from dagmcwritemcnp_
/////////////
//////////////////////////////////////////////////////////////////////////
// void dagmcwritemcnp_(char *lfile, int *llen)  // file with cell/surface cards
void fludagwrite(std::string fname )  // file with cell/surface cards
{
  MBErrorCode rval;

  // lfile[*llen]  = '\0';

  std::vector< std::string > fluka_keywords;
  std::map< std::string, std::string > fluka_keyword_synonyms;

  fluka_keywords.push_back( "mat" );
  fluka_keywords.push_back( "rho" );
  fluka_keywords.push_back( "comp" );
  fluka_keywords.push_back( "imp.n" );
  fluka_keywords.push_back( "imp.p" );
  fluka_keywords.push_back( "imp.e" );
  fluka_keywords.push_back( "tally" );
  fluka_keywords.push_back( "spec.reflect" );
  fluka_keywords.push_back( "white.reflect" );
  fluka_keywords.push_back( "graveyard" );
  
  // fluka_keyword_synonyms[ "rest.of.world" ] = "graveyard";
  // fluka_keyword_synonyms[ "outside.world" ] = "graveyard";

  // parse data from geometry
  rval = DAG->parse_properties (fluka_keywords);
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  std::cerr << "Going to write an lcad file = " << fname << std::endl;
  // Before opening file for writing, check for an existing file
  if( fname != "lcad" ){
    // Do not overwrite a lcad file if it already exists, except if it has the default name "lcad"
    if( access( fname.c_str(), R_OK ) == 0 ){
      std::cout << "DagMC: reading from existing lcad file " << fname << std::endl;
      return; 
    }
  }

  std::ofstream lcadfile( fname.c_str() );

  int num_cells = DAG->num_entities( 3 );
  int num_surfs = DAG->num_entities( 2 );

  int cmat = 0;
  double crho, cimp = 1.0;

  // write the cell cards
  for( int i = 1; i <= num_cells; ++i ){

    MBEntityHandle vol = DAG->entity_by_index( 3, i );
    int cellid = DAG->id_by_index( 3, i );
    // set default importances for p and e to negative, indicating no card should be printed.
    double imp_n = 1, imp_p = -1, imp_e = -1;

    if( DAG->has_prop( vol, "imp.n" )){
      get_real_prop( vol, cellid, "imp.n", imp_n );
    }

    if( DAG->has_prop( vol, "imp.p" )){
      get_real_prop( vol, cellid, "imp.p", imp_p );
    }

    if( DAG->has_prop( vol, "imp.e" )){
      get_real_prop( vol, cellid, "imp.e", imp_e );
    }

    lcadfile << cellid << " ";

    bool graveyard = DAG->has_prop( vol, "graveyard" );

    if( graveyard ){
      lcadfile << " 0 imp:n=0";
      if( DAG->has_prop(vol, "comp") ){
        // material for the implicit complement has been specified.
        get_int_prop( vol, cellid, "mat", cmat );
        get_real_prop( vol, cellid, "rho", crho );
        std::cout << "Detected material and density specified for implicit complement: " << cmat <<", " << crho << std::endl;
        cimp = imp_n;
      }
    }
    else if( DAG->is_implicit_complement(vol) ){
      lcadfile << cmat;
      if( cmat != 0 ) lcadfile << " " << crho;
      lcadfile << " imp:n=" << cimp << " $ implicit complement";
    }
    else{
      int mat = 0;
      get_int_prop( vol, cellid, "mat", mat );

      if( mat == 0 ){
        lcadfile << "0";
      }
      else{
        double rho = 1.0;
        get_real_prop( vol, cellid, "rho", rho );
        lcadfile << mat << " " << rho;
      }
      lcadfile << " imp:n=" << imp_n;
      if( imp_p > 0 ) lcadfile << " imp:p=" << imp_p;
      if( imp_e > 0 ) lcadfile << " imp:e=" << imp_e;
    }

    lcadfile << std::endl;
  } // end iteration through cells

  // cells finished, skip a line
  lcadfile << std::endl;
  
  // jcz - do we need this?
  // write the surface cards
/*
  for( int i = 1; i <= num_surfs; ++i ){
    MBEntityHandle surf = DAG->entity_by_index( 2, i );
    int surfid = DAG->id_by_index( 2, i );

    if( DAG->has_prop( surf, "spec.reflect" ) ){
      lcadfile << "*";
    }
    else if ( DAG->has_prop( surf, "white.reflect" ) ){
      lcadfile << "+";
    }
    lcadfile << surfid << std::endl;
  }

  // surfaces finished, skip a line
  lcadfile << std::endl;
*/
  // write the tally cards
  std::vector<std::string> tally_specifiers;
  rval = DAG->get_all_prop_values( "tally", tally_specifiers );
  if( rval != MB_SUCCESS ) exit(EXIT_FAILURE);

  for( std::vector<std::string>::iterator i = tally_specifiers.begin();
       i != tally_specifiers.end(); ++i )
  {
    int dim = 0;
    char* card = get_tallyspec( *i, dim );
    if( card == NULL ){
      std::cerr << "Invalid dag-mcnp tally specifier: " << *i << std::endl;
      std::cerr << "This tally will not appear in the problem." << std::endl;
      continue;
    }
    std::stringstream tally_card;

    tally_card << card;
    std::vector<MBEntityHandle> handles;
    std::string s = *i;
    rval = DAG->entities_by_property( "tally", handles, dim, &s );
    if( rval != MB_SUCCESS ) exit (EXIT_FAILURE);

    for( std::vector<MBEntityHandle>::iterator j = handles.begin();
         j != handles.end(); ++j )
    {
      tally_card << " " << DAG->get_entity_id(*j);
    }

    tally_card  << " T";
    delete[] card;

    // write the contents of the the tally_card without exceeding 80 chars
    std::string cardstr = tally_card.str();
    while( cardstr.length() > 72 ){
        size_t pos = cardstr.rfind(' ',72);
        lcadfile << cardstr.substr(0,pos) << " &" << std::endl;
        lcadfile << "     ";
        cardstr.erase(0,pos);
    }
    lcadfile << cardstr << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////
///////////////
//////////////
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

