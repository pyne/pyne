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
#include "UnitNumberManager.hpp"

#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "DagMC.hpp"
#include "moab/Types.hpp"
using moab::DagMC;

#include <iomanip>
#include <sstream>
#include <set>
#include <cstring>
#include <string>

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

std::multimap<std::string, unsigned int> scoring_vol_map;

UnitNumberManager unit_no_mgr = UnitNumberManager();

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

//----------------------------------*-C++-*----------------------------------//
/* Called by mainReadVol
 *  cpp_dagmcinit is called directly from c++ or from a fortran-called wrapper.
 *  Precondition:  myfile exists and is readable
 * \file   ~/DAGMC/FluDAG/src/cpp/mainReadVol.cpp
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
// - curRegion is not currently used.  
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
///////			End normal() and f_normal()
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
// Question:  Should pV, the direction vector, be used?  
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

void lkmgwr(double& pSx, double& pSy, double& pSz,
            double* pV, const int& oldReg, const int& oldLttc,
	    int& flagErr, int& newReg, int& newLttc)
{
    std::cerr<<"============= LKMGWR =============="<<std::endl;
    const double xyz[] = {pSx, pSy, pSz}; // location of the particle (xyz)
    int is_inside = 0; // logical inside or outside of volume
    int num_vols = DAG->num_entities(3); // number of volumes

    for (int i = 1 ; i <= num_vols ; i++) // loop over all volumes
      {
	MBEntityHandle volume = DAG->entity_by_index(3, i); // get the volume by index
	// No ray history or ray direction.
	MBErrorCode code = DAG->point_in_volume(volume, xyz, is_inside);

	// check for non error
	if(MB_SUCCESS != code) 
	  {
	    std::cout << "Error return from point_in_volume!" << std::endl;
	    flagErr = 1;
	    return;
	  }

	if ( is_inside == 1 ) // we are inside the cell tested
	  {
	    newReg = i;
	    flagErr = i+1;
	    if(debug)
	      {
		std::cout << "point is in region = " << newReg << std::endl;
	      }
	    return;
	  }
      }  // end loop over all volumes

    std::cout << "particle is nowhere!" << std::endl;
    newReg = -100;
    std::cout << "point is not in any volume" << std::endl;
    return;
}


//--------------------------------------------------------------------------//
// f_lookdb(..)
//---------------------------------------------------------------------------//
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

//--------------------------------------------------------------------------//
// f_g1rt()
//---------------------------------------------------------------------------//
void f_g1rt(void)
{
  if(debug)
    {
      std::cout<<"============ F_G1RT ============="<<std::endl;
    }
    return;
}

//---------------------------------------------------------------------------//
// f_idnr(..)
//---------------------------------------------------------------------------//
// Set DNEAR option if needed
int f_idnr(const int & nreg, const int & mlat) 
{
   // returns 0 if user doesn't want Fluka to use DNEAR to compute the 
   // step (the same effect is obtained with the GLOBAL (WHAT(3)=-1)
   // card in fluka input), returns 1 if user wants Fluka always to use DNEAR.

   return 0;
}

//---------------------------------------------------------------------------//
// rg2nwr(..)
//---------------------------------------------------------------------------//
// Wrapper for getting region name corresponding to given region number
void rg2nwr(const int& mreg, const char* Vname)
{
  std::cerr << "============= RG2NWR ==============" << std::endl;    
  std::cerr << "mreg=" << mreg << std::endl;
  char * vvname;
  region2name(mreg, vvname);
  Vname = vvname;
  std::cerr << "reg2nmwr: Vname " << Vname<< std::endl;  
  return;
}

/**************************************************************************************************/
/******                                End of FLUKA stubs                                  ********/
/**************************************************************************************************/

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
//---------------------------------------------------------------------------//
// fludagwrite_assignma
//---------------------------------------------------------------------------//
/// Called from mainFludag when only one argument is given to the program.
//  This function writes out a simple numerical material assignment to the named argument file
//  Example usage:  mainFludag dagmc.html
//  Outputs
//           mat.inp  contains MATERIAL and ASSIGNMAt records for the input geometry.
//                    The MATERIAL is gotten by parsing the Cubit volume name on underscores.  
//                    The string after "M_" is considered to be the material for that volume.
//                    There are no MATERIAL cards for the materials in the FLUKA_mat_set list
//                    For the remaining materials, there is one MATERIAL card apiece (no dups)
//                    User-named (not predefined) materials are TRUNCATED to 8 chars.
//                    User-named material id's start at 25 and increment by 1 for each MATERIAL card
//           index-id.txt  Map of FluDAG volume index vs Cubit volume ids, for info only.
//  Note that a preprocessing step to this call sets up the the DAG object that contains 
//  all the geometry information contained in dagmc.html.  
//  the name of the (currently hardcoded) output file is "mat.inp"
//  The graveyard is assumed to be the last region.
void fludagwrite_assignma(std::string filename_to_write)  // file with cell/surface cards
{
  int num_vols = DAG->num_entities(3);
  // std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
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

  // jcz debug: lists DEN, M, NEUTRON, S, USRTRACK
  std::vector< std::string >::iterator vit;
  for (vit=keywords.begin(); vit!=keywords.end(); ++vit)
  {
    std::cout << *vit << std::endl;
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
  std::ostringstream graveyard_str;
  std::ostringstream impl_compl_str;

  // open an outputstring for the M_  (ASSIGNMAt) portions
  std::ostringstream A_filestr;

  // Open an outputstring for index-id table and put a header in it
  std::ostringstream idstr;
  idstr << std::setw(5) <<  "Index" ;
  idstr << std::setw(5) <<  "   Id" << std::endl;

  // Prepare a list to contain unique materials not in Fluka's list
  std::list<std::string> uniqueMatList;

  // Loop through 3d entities (vols).  
  std::vector<std::string> vals;
  std::string material_trunc;
  char buffer[MAX_MATERIAL_NAME_SIZE];
  for (unsigned int i=1; i<=num_vols ; i++)
  {  
      vals.clear();
      entity = DAG->entity_by_index(3, i);

      // Create the id-index string for this vol
      addToIDIndexFile(i, idstr);

      // Create the mat.inp string for this vol
      if (DAG->has_prop(entity, "graveyard"))
      {
	 graveyard_str << std::setw(10) << std::left  << "ASSIGNMAt";
	 graveyard_str << std::setw(10) << std::right << "BLCKHOLE";
	 graveyard_str << std::setw(10) << std::right << i << std::endl;
      }
      if (DAG->has_prop(entity, "M"))
      {
         DAG->prop_values(entity, "M", vals);

         process_Mi(A_filestr, entity, uniqueMatList, i);
      } // end processing of "M_" property
      if (DAG->has_prop(entity, "S"))
      {
         process_Si(entity, i);
      } // end processing of "S_" property
  }  // end of volume processing loop

  std::ostringstream S_filestr;
  std::ostringstream r_filestr;
  std::ostringstream ut_filestr;
  std::ostringstream uc_filestr;
  std::ostringstream ub_filestr;
  std::ostringstream uy_filestr;

  // Print out the scoring by volume map
  std::cout << "All scoring.particle and volumes" << std::endl;
  unsigned int counter = 0;
  std::map<std::string, unsigned int>::iterator uit;
  
  // Go through the map that was created while going through the volumes
  for (uit = scoring_vol_map.begin(); uit != scoring_vol_map.end(); ++uit)
  {
     counter++;
     // std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
     std::cout << counter << ". " << uit->first << " => " << uit->second << std::endl;

     // Get the volume id of the current volume, whose scoring info we are pulling out
     int iVol = uit->second;

     // Create a vector of '.'-delimited strings from the original "key" part of the 
     // group name (now the "value" part of the map).  
     // The vector may be of size 1, 2, or 3, depending on the score type;
     // No more than the first two values are used in the unit map.
     std::vector<std::string> score_words = StringSplit(uit->first,".");
     std::string score_name;
     char particle[3];
     
     if (score_words.size() == 0)
     {
        std::cout << "The score group for volume " << iVol << " is empty." << std::endl;
        continue;
     }
     if (score_words.size() > 0)
     {
        score_name = score_words[0];
     }
     if (score_words.size() > 1)
     {

        std::size_t len = score_words[1].copy(particle,2);
        particle[len] = '\0';
     }
    
     char strDetName[10];
     float fortran_unit;
     double measurement;

     // The RESNUCLEI section:  use to_string when c11 comes
     if (score_name.compare("RESNUCLEI") == 0)
     {
        if (r_filestr.tellp() == 0)
        {
           r_filestr << "* RESNUCLEI scoring requests." << std::endl;
           r_filestr << header << std::endl;
        }
        // jcz ToDo;  Error return.  How serious an error is -1?  Is returning warranted?
        measurement = measurementOfVol(iVol);
        fortran_unit = get_score_particle_unit(score_words);
        sprintf (strDetName, "%s%d","RES_", iVol);

        r_filestr << std::setw(10) << std::left << "RESNUCLEI";
        r_filestr << std::setw(20) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
	r_filestr << std::setw(30) << std::right << std::fixed  << std::setprecision(1) << (float)iVol;
        r_filestr << std::setw(9)  << std::right << measurement << " ";
        r_filestr << std::setw(10) << std::left <<  strDetName;

        r_filestr << std::endl;
     }
     // The USRTRACK and USRCOLL section 
     else if (score_name.compare("USRTRACK") == 0 || score_name.compare("USRCOLL") == 0)
     {
        // These tallies need both a score and a particle name
	if (score_words.size() >= 2)
        {
           fortran_unit = get_score_particle_unit(score_words);
           if (score_name.compare("USRTRACK") == 0)
           {
              if (ut_filestr.tellp() == 0)
              {
                 ut_filestr << "* USRTRACK scoring requests." << std::endl;
                 ut_filestr << header << std::endl;
              }
              sprintf (strDetName, "%s_%s_%d","TR",  particle, iVol);
              basic_score(ut_filestr, score_words, iVol, fortran_unit, measurement, strDetName);
           }
 	   else  // USRCOLL case
           {
              if (uc_filestr.tellp() == 0)
              {
                 uc_filestr << "* USRCOLL scoring requests." << std::endl;
                 uc_filestr << header << std::endl;
              }
              sprintf (strDetName, "%s_%s_%d","CO",  particle, iVol);
              basic_score(uc_filestr, score_words, iVol, fortran_unit, measurement, strDetName);
           }
        }
        else   // Error:  no particle was added to the group name
        {
           std::cerr << "Error: the " << score_name << " score does not include a particle. " 
                     << "Please label the group with particle" << std::endl;
        }
     } 
     // The USRBDX section
     else if (score_name.compare("USRBDX") == 0 || score_name.compare("USRYIELD") == 0)
     {
        // These tallies need a score, particle name, and a FROM volume
        if (score_words.size() >= 3)
        {
           fortran_unit = get_score_particle_unit(score_words);
           if (score_name.compare("USRBDX") == 0)     // USRBDX case 
           {
              if (ub_filestr.tellp() == 0)
              {
                 ub_filestr << "* USRBDX scoring requests." << std::endl;
                 ub_filestr << header << std::endl;
              }
              two_vol_score(ub_filestr, score_words, iVol, fortran_unit, "BDX");
           }
           else   // USRYIELD case
           {
              if (uy_filestr.tellp() == 0)
              {
                 uy_filestr << "* USRYIELD scoring requests." << std::endl;
                 uy_filestr << header << std::endl;
              }
              two_vol_score(uy_filestr, score_words, iVol, fortran_unit, "YIELD");
           }
        }
        else   // Error:  missing particle or volume
        {
           std::cerr << "Error: the " << score_name << " score is missing a particle or a volume. " 
                     << "Please label the group with both particle and FROM volume" << std::endl;
        }
     } 
  }

  // Collect all the scoring records into one stream
  S_filestr <<  r_filestr.str();
  S_filestr << ut_filestr.str();
  S_filestr << uc_filestr.str();
  S_filestr << ub_filestr.str();
  S_filestr << uy_filestr.str();
  // Optional: send all the scores to the screen
  std::cout <<  S_filestr.str();

  // Finish the ostr with the implicit complement card
  std::string implicit_comp_comment = "* The next volume is the implicit complement";
  impl_compl_str << implicit_comp_comment << std::endl;
  impl_compl_str << std::setw(10) << std::left  << "ASSIGNMAt";
  impl_compl_str << std::setw(10) << std::right << "VACUUM";
  impl_compl_str << std::setw(10) << std::right << num_vols << std::endl;

  // Prepare the MATERIAL cards as a string stream using a list
  // of materials that has no duplicates
  uniqueMatList.sort();
  uniqueMatList.unique();
  std::ostringstream MAT_filestr;
  processUniqueMaterials(MAT_filestr, uniqueMatList, header);
  // Print the final list
  if (debug)
  {
     std::list<std::string>::iterator it; 
     for (it=uniqueMatList.begin(); it!=uniqueMatList.end(); ++it)
     {
        std::cout << *it << std::endl;
     }
  
     // Show the output string just created
     std::cout << ostr.str();
  }

  // Prepare a file of the given name; put a header in it;  This file contains the FLUKA input cards 
  // (records) the user has requested via geometry tags.
  std::ofstream records_filestr( filename_to_write.c_str());
  records_filestr << header << std::endl;

  // Put all the filestr parts together
  records_filestr << MAT_filestr.str();    // The material list is created separately
  records_filestr << graveyard_str.str();  // the graveyard
  records_filestr << A_filestr.str();      // ASIGNMAt statements
  records_filestr << impl_compl_str.str(); // implicit complement
  records_filestr << S_filestr.str();      // Detector tallies
  records_filestr.close();

  std::cout << "Writing input file = " << filename_to_write << std::endl; 

  // Prepare an output file named "index_id.txt" for idstr
   writeToFileNamed(idstr, "index_id.txt");
}
// End fludagwrite_assignma

//---------------------------------------------------------------------------//
// measurementOfVol
//---------------------------------------------------------------------------//
// Convenience function to call dag measure_volume(..)
double measurementOfVol(int iVol)
{
   MBErrorCode ret;

   // Calculate some things we'll need that are based on the volume id 
   EntityHandle handle = DAG->entity_by_id(3, iVol);

   double measurement;
   ret = DAG->measure_volume(handle, measurement);
   if (MB_SUCCESS != ret) 
   {
      std::cerr << "DAGMC failed to get the measured volume of region " <<  iVol <<  std::endl;
      measurement = -1.0;
   }
   return measurement; 
}
//---------------------------------------------------------------------------//
// usrtrackRecord
//---------------------------------------------------------------------------//
// Specialized record-preparer for USRTRACK defaults
/*
void usrtrackRecord(std::ostringstream& ostr, 
                 std::vector<std::string> score_words, 
                 int iVol, float fortran_unit, std::string score_prefix)
{
     // Prepare the detector name to go at the end of the line
     char strDetName[10];
     std::string subname (score_prefix + "_" + score_words[2] + "_");
     char *cstr = new char [subname.length() + 1];
     std::strcpy (cstr, subname.c_str());
     sprintf (strDetName, "%s%d", cstr, iVol);

     // We are guaranteed there are three values in score_words
     ostr << std::setw(10) << std::left << score_words[0];           
     ostr << std::setw(10) << std::right << "-1.0";
     ostr << std::setw(10) << std::right << score_words[1];           
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << score_words[2];
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << (float)iVol;

     ostr << std::setw(10) << std::right << " ";
     ostr << std::setw(10) << std::left <<  strDetName << std::endl;

     sdum_endline(ostr, score_words[0]);
     if (score_words.size() > 3)  // hmm, there is a particle piece, but also more
     {
           std::cerr << "Error:  the " << score_words[0] << " score has more than one particle reference.  " 
                     << "Only the first particle, " << score_words[1] << ", is used." << std::endl;  
     }
}
*/
//---------------------------------------------------------------------------//
// two_vol_score
//---------------------------------------------------------------------------//
// Some records are very similar, differing only in name
void two_vol_score(std::ostringstream& ostr, 
                 std::vector<std::string> score_words, 
                 int iVol, float fortran_unit, std::string score_prefix)
{
     // Prepare the detector name to go at the end of the line
     char strDetName[10];
     std::string subname (score_prefix + "_" + score_words[2] + "_");
     char *cstr = new char [subname.length() + 1];
     std::strcpy (cstr, subname.c_str());
     sprintf (strDetName, "%s%d", cstr, iVol);

     // We are guaranteed there are three values in score_words
     ostr << std::setw(10) << std::left << score_words[0];           
     ostr << std::setw(10) << std::right << "-1.0";
     ostr << std::setw(10) << std::right << score_words[1];           
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << score_words[2];
     ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << (float)iVol;

     ostr << std::setw(10) << std::right << " ";
     ostr << std::setw(10) << std::left <<  strDetName << std::endl;

     sdum_endline(ostr, score_words[0]);
     if (score_words.size() > 3)  // hmm, there is a particle piece, but also more
     {
           std::cerr << "Error:  the " << score_words[0] << " score has more than one particle reference.  " 
                     << "Only the first particle, " << score_words[1] << ", is used." << std::endl;  
     }
}
//---------------------------------------------------------------------------//
// basic_score
//---------------------------------------------------------------------------//
// Some records are very similar, differing only in name
// This handls USRCOLL and USRTRACK score requests
void basic_score(std::ostringstream& ostr, 
                 std::vector <std::string> score_words, 
                 int iVol, float fortran_unit, double measurement,
                 std::string name)
{
     measurement = measurementOfVol(iVol);
     if (score_words.size() >= 2)  // all is good
     {
           ostr << std::setw(10) << std::left << score_words[0];           
           ostr << std::setw(10) << std::right << "-1.0";
           ostr << std::setw(10) << std::right << score_words[1];           
           ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << fortran_unit;
	   ostr << std::setw(10) << std::right << std::fixed << std::setprecision(1) << (float)iVol;
           ostr << std::setw(10) << std::right << measurement << " ";
           // Default number of energy bins
           ostr << std::setw(8)  << std::right << "1.";

           // ostr << std::setw(10) << std::left <<  name << std::endl;;
           ostr << " " << name << std::endl;;
	   sdum_endline(ostr, score_words[0], true);
     }
     if (score_words.size() > 2)  // hmm, there is a particle word, but also more
     {
           std::cerr << "Warning:  the " << score_words[0] << " score has more than one particle reference.  " 
                     << "Only the first particle, " << score_words[1] << ", is used." << std::endl;  
     }
}
//---------------------------------------------------------------------------//
// sdum_endline
//---------------------------------------------------------------------------//
// Create a standard continuation line, putting the '&' on the 71st space
// The line is tacked on to the output stream
void sdum_endline(std::ostringstream& ostr, std::string score, bool isBasic)
{
     ostr << std::setw(10) << std::left << score;           
     if (isBasic)  // This path is for using the defaults common to USRTRACK and USRYIELD
     {
        // Default maximum energy
        ostr << std::setw(10) << std::right << "1E20";           
        // Default minimum energy
        ostr << std::setw(10) << std::right << "1.0";           
        ostr << std::setw(41) << std::right << "&";
     }
     else
     {
        ostr << std::setw(61) << std::right << "&";
     }
     ostr << std::endl;
}

//---------------------------------------------------------------------------//
// process_Si
//---------------------------------------------------------------------------//
// Process group names that have S as the key and track.[p'le] as the value
// Examples:  USRTRACK.PROTON, RESNUCLEI, USRCOLL.NEUTRON
// A multimap is used to store every score request for every volume,
void process_Si(MBEntityHandle entity, unsigned int vol_id)
{
    MBErrorCode ret;
    std::vector<std::string> vals;

    // We only get here if has_prop(... "S" ...) is true
    ret = DAG->prop_values(entity, "S", vals);
    if (MB_SUCCESS != ret) 
    {
       std::cerr << "DAGMC failed to get S_ properties" <<  std::endl;
       return;
    }
    for (int i=0; i<vals.size(); i++) 
    { 
        scoring_vol_map.insert(std::pair<std::string, unsigned int>(vals[i], vol_id)); 
    }
}

//---------------------------------------------------------------------------//
// get_score_particle_mapname
//---------------------------------------------------------------------------//
/// Convert the first two values of the group name (in the case there are two) to a 
//  standard name for the fortran unit number getter
std::string get_score_particle_mapname(std::string score_name, std::string particle_name)
{
    std::string keyword (score_name + "." + particle_name);
    char *cstr = new char [keyword.length() + 1];
    std::strcpy (cstr, keyword.c_str());
    return cstr;
}

//---------------------------------------------------------------------------//
// get_score_particle_unit
//---------------------------------------------------------------------------//
/// Parse the vector of strings from the scoring group name values and determine
//  the correct fortran unit number for them.  
//  There can be 1, 2, or more words in teh score_words vectore. 
//  If 1, it is used.  If 2, both are used.  If more than 2, only the first two are used.
// This function relies on the UnitNumberManager class, whose key method returns
// an int, however we always need a float unit number for the fluka cards, so
// a float is returned. 
//  If the vector is empty, this function returns -1
float get_score_particle_unit(std::vector<std::string> score_words)
{
    std::string mapname;
    if (score_words.size() == 1)
    {
        mapname = score_words[0];
    }
    else if (score_words.size() > 1)
    {
        mapname = get_score_particle_mapname(score_words[0], score_words[1]);
    }
    return (float)unit_no_mgr.getUnitNumber(mapname);
}
//---------------------------------------------------------------------------//
// processUniqueMaterials
//---------------------------------------------------------------------------//
// Convenience method to create MATERIAL cards
void processUniqueMaterials(std::ostringstream& ostr, 
                            std::list<std::string> uniqueList, 
                            std::string header)
{
  // Prepare the MATERIAL cards as a string stream
  if (uniqueList.size() != 0)
  {
     int matID = 25;
     std::list<std::string>::iterator it; 
     for (it=uniqueList.begin(); it!=uniqueList.end(); ++it)
     {
        ostr << std::setw(10) << std::left << "MATERIAL";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << ++matID;
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::right << "";
        ostr << std::setw(10) << std::left << *it << std::endl;
     }
  }
  if (uniqueList.size() !=0)
  {
     ostr << header << std::endl;
  }
  return;
}
//---------------------------------------------------------------------------//
// process_MI
//---------------------------------------------------------------------------//
// Add a record to ostr of the form "ASSIGNMA ....."
void process_Mi(std::ostringstream& ostr, 
                MBEntityHandle entity, 
                std::list<std::string> &matList, unsigned i)
{
    MBErrorCode ret;
    std::vector<std::string> vals;
    char buffer[MAX_MATERIAL_NAME_SIZE];
    std::string material_trunc;

    ret = DAG->prop_values(entity, "M", vals);
    if (MB_SUCCESS != ret) 
    {
       std::cerr << "DAGMC failed to get M_ properties" <<  std::endl;
       return;
    }

    if (vals.size() >= 1)
    {
       // Make a copy of string in vals[0]; full string needs to be compared to
       // FLUKA materials list; copy is for potential truncation
       std::strcpy(buffer, vals[0].c_str());
       material_trunc = std::string(buffer);
            
       if (vals[0].size() > 8)
       {
           material_trunc.resize(8);
       }

       if (FLUKA_mat_set.find(vals[0]) == FLUKA_mat_set.end())
       {
          // current material is not in the pre-existing FLUKA material list
          matList.push_back(material_trunc); 
          std::cout << "Adding material " << material_trunc << " to the MATERIAL card list" << std::endl;
       }
    }
    else
    {
         material_trunc = "moreThanOne";
    }
     ostr << std::setw(10) << std::left  << "ASSIGNMAt";
     ostr << std::setw(10) << std::right << material_trunc;
     ostr << std::setw(10) << std::right << i << std::endl;
}
//---------------------------------------------------------------------------//
// addToIDIndexMap(int i, 
//---------------------------------------------------------------------------//
// Convenience method to connect the geometry id to the ith volume 
void addToIDIndexFile(int i, std::ostringstream &idstr)
{
      idstr << std::setw(5) << std::right << i;
      idstr << std::setw(5) << std::right << DAG->id_by_index(3,i) << std::endl;
}
//---------------------------------------------------------------------------//
// writeStringToFile
//---------------------------------------------------------------------------//
// Convenience method to write a prepared stream to a file
void writeToFileNamed(std::ostringstream& file_contents, std::string filename)
{
     std::ofstream file_stream(filename.c_str());
     file_stream << file_contents.str();
     file_stream.close(); 
}

//---------------------------------------------------------------------------//
// process_SI
//---------------------------------------------------------------------------//
// Process the request for a RESNUCLEI tally
/*
void process_S_RESNUCLEI(std::ostringstream& ostr, MBEntityHandle entity, unsigned i)
{
   MBErrorCode ret;
   std::vector<std::string> vals;

   int unitNum = getRESNUCLEIUnitNumber()
   return;
} 
*/
/*
int getRESNUCLEIUnitNumber()
{
    if (RESNUCLEI_Unit = -1)
    {
       RESNUCLEI_Unit = getNextUnitNumber();
       ++num_units_in_use;
    }
    return RESNUCLEI_Unit;
}
*/
//---------------------------------------------------------------------------//
// process_MI
//---------------------------------------------------------------------------//
// 
void process_Mi(std::ostringstream& ostr, MBEntityHandle entity, std::list<std::string> &matList, unsigned i)
{
    MBErrorCode ret;
    std::vector<std::string> vals;
    char buffer[MAX_MATERIAL_NAME_SIZE];
    std::string material_trunc;

    ret = DAG->prop_values(entity, "M", vals);
    if (MB_SUCCESS != ret) 
    {
       std::cerr << "DAGMC failed to get M_ properties" <<  std::endl;
       return;
    }

    if (vals.size() >= 1)
    {
       // Make a copy of string in vals[0]; full string needs to be compared to
       // FLUKA materials list; copy is for potential truncation
       std::strcpy(buffer, vals[0].c_str());
       material_trunc = std::string(buffer);
            
       if (vals[0].size() > 8)
       {
           material_trunc.resize(8);
       }

       if (FLUKA_mat_set.find(vals[0]) == FLUKA_mat_set.end())
       {
          // current material is not in the pre-existing FLUKA material list
          matList.push_back(material_trunc); 
          std::cout << "Adding material " << material_trunc << " to the MATERIAL card list" << std::endl;
       }
     }
     else
     {
         material_trunc = "moreThanOne";
     }
     ostr << std::setw(10) << std::left  << "ASSIGNMAt";
     ostr << std::setw(10) << std::right << material_trunc;
     ostr << std::setw(10) << std::right << i << std::endl;
}
//---------------------------------------------------------------------------//
// getNextUnitNumber()
//---------------------------------------------------------------------------//
// Convenience method to get the next logical unit number for the writing-out 
// field of a FLUKA card.  The key is when to call.  
/*
int getNextUnitNumber()
{
    int retval =  START_UNIT - num_units_in_use;
    ++num_units_in_use;
    return retval;
}
*/
//---------------------------------------------------------------------------//
// addToIDIndexMap(int i, 
//---------------------------------------------------------------------------//
// Convenience method to connect the geometry id to the ith volume 
void addToIDIndexFile(int i, std::ostringstream &idstr)
{
      idstr << std::setw(5) << std::right << i;
      idstr << std::setw(5) << std::right << DAG->id_by_index(3,i) << std::endl;
}
//---------------------------------------------------------------------------//
// writeStringToFile
//---------------------------------------------------------------------------//
// Convenience method to write a prepared stream to a file
void writeToFileNamed(std::ostringstream& file_contents, std::string filename)
{
     std::ofstream file_stream(filename.c_str());
     file_stream << file_contents.str();
     file_stream.close(); 
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
//                      Called in WrapReg2Name
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
void dagmc_version_(double* dagmcVersion)
{
  *dagmcVersion = DAG->version();
}

