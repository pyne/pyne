#include "fluka_funcs.h"

#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "DagMC.hpp"
using moab::DagMC;

#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

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

std::string ExePath() {
    int MAX_PATH = 256;
    char buffer[MAX_PATH];
    getcwd(  buffer, MAX_PATH );
    // std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
    // return std::string( buffer ).substr( 0, pos);
    return std::string( buffer );
}

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
  std::string str1="../";
  std::string str2= std::string(cfile);
  str1.append(str2);
  std::cout << "\nmy file is " << str1 << "\n" << std::endl;
  char *myfile;
  myfile = &str1[0];
  

    // initialize this as -1 so that DAGMC internal defaults are preserved
    // user doesn't set this
  double arg_facet_tolerance = -1;
                                                                        
  // jcz: leave arg_facet_tolerance as defined on previous line
  // if ( *ftlen > 0 ) arg_facet_tolerance = atof(ftol);

  // read geometry
  // rval = DAG->load_file(cfile, arg_facet_tolerance );
  rval = DAG->load_file(myfile, arg_facet_tolerance );
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


/**************************************************************************************************/
/******                                FLUKA stubs                                         ********/
/**************************************************************************************************/
extern "C" int lookdb_ (double *X, double *Y, double *Z, int *numErrLm)
{
	std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << std::endl;
	return *numErrLm;
}
extern "C" int lookmg_ (double *X, double *Y, double *Z,
                  double dir[3], // Direction cosines vector
		  int *RegionNum, // region number
                  int *newCell,    // Output: region # of p'le after step ("IRPRIM" in FLUKA)
                  int *Ierr        // Output: error code
)
{
	std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << std::endl;
	return *Ierr;
}

extern "C" int lookfx_ (double *X, double *Y, double *Z,
                  double dir[3], // Direction cosines vector
		  int *RegionNum, // region number
                  int *newCell,    // Output: region # of p'le after step ("IRPRIM" in FLUKA)
                  int *Ierr        // Output: error code
)
{
	std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << std::endl;
	return *Ierr;
}

/*
 * Allows tracking initialization.
 * Same input and output as LOOKMG.
*/
extern "C" int lookz_ (double *X, double *Y, double *Z,
                  double dir[3], // Direction cosines vector
		  int *RegionNum, // region number
                  int *newCell,    // region # of p'le after step ("IRPRIM" in FLUKA)
                  int *Ierr        // error code
)
{
	std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << std::endl;
        std::cout << "In C++ function lookz_" << std::endl;
	std::cout << "\tINPUT: Position X,Y,Z  " << *X << ", " \
                           <<  *Y  << ", "  << *Z << std::endl;
	std::cout << "\tINPUT: Direction cosines   " << dir[1] << ", " \
                           <<  dir[2]  << ", "  << dir[3] << std::endl;
	std::cout << "\tINPUT: RegionNum " << *RegionNum << std::endl;


        *newCell = *RegionNum + 1;
	std::cout << "OUTPUT: int *newCell = " << *newCell << std::endl;
	return *Ierr;
}

/*
 * Returns a unit vector at the previous point of the tracking (in case a boundary
 * crossing occurred), that is the intersection point between the path and the
 * preceding boundary.  The call is done from the routine GEONOR and all
 * cosines must be normalized at 10e-16, and they must be consistent with the
 * actual normal, with similar accuracy.
*/
extern "C" int norml_(double *U, double *V, double *W)
{
	
	std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << std::endl;
	std::cout << " Norml unit vector: U, V, W = " << *U << ", " \
                           <<  *V  << ", "  << *W << std::endl;
        *U = 98.7;
	*V = 6e-5;
	*W = 4.3e2;	
	return 1;
}
/**
*
*  From "Building an Interface between FLUKA and the GEANT4 Geometry Package":
*  G1FLU calculates the distance travelled in the present zone/region and the 
*  number of the next zone/region to be entered by the particle.
*  NOTE:  All the variables are in double precision.
* 
*  Input
*   int *SurfaceNum Surface number hit by the p'le at the preceding step 
*                   (-1 if the p'le changed its direction or if it's a new
*                    particle)
* Output
*  int *SurfaceNum  number of surface hit by the particle (1 for normal
*                   tracking, 0 otherwise      
*  double *tol      parameter with value less than or equal to the 
*                   MINIMUM of the distances between the particle and each
*                   boundary (if the boundary surface is too complex 0 is returned)
*
* Return
*  double *tol      Per the documentation G1FLU returns DSNEAR
*/
extern "C" double g1flu_(double pos[3], // Cartesian coordinate vector
                  double dir[3], // Direction cosines vector
		  int *RegionNum, // region number
                  double *dist,   // step length SUGGESTED
                  int *SurfaceNum, // "NASC" in FLUKA, see above for text
                  double *step,    // step length approved
                  int *newCell,    // region # of p'le after step ("IRPRIM" in FLUKA)
                  double *tol,     // "DSNEAR" in FLUKA, see above for text
                  int *Ierr        // error code
)
{
	std::cout << __FILE__ << ", " << __func__ << ":" << __LINE__ << std::endl;
	std::cout << " Cartesian coordinate vector: pos[3] = " << pos[0]<< ", " \
                           <<  pos[1]  << ", "  << pos[2] << std::endl;
	double huge = 10e10;
	double dls;
        int jap;
	int jsu;
	int nps;

/* ToDo:  
        dagmctrack_(RegionNum, &(dir[0]), &(dir[1]), &(dir[2]), &(pos[0]), &(pos[1]), &(pos[2]), 
                     &huge, &dls, &jap, &jsu, &nps);
*/
	return *tol;
}

/**************************************************************************************************/
/******                                End of FLUKA stubs                                  ********/
/**************************************************************************************************/

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

void dagmcwritemcnp_(char *lfile, int *llen)  // file with cell/surface cards
                     
{
  MBErrorCode rval;

  lfile[*llen]  = '\0';

  std::vector< std::string > mcnp5_keywords;
  std::map< std::string, std::string > mcnp5_keyword_synonyms;

  mcnp5_keywords.push_back( "mat" );
  mcnp5_keywords.push_back( "rho" );
  mcnp5_keywords.push_back( "comp" );
  mcnp5_keywords.push_back( "imp.n" );
  mcnp5_keywords.push_back( "imp.p" );
  mcnp5_keywords.push_back( "imp.e" );
  mcnp5_keywords.push_back( "tally" );
  mcnp5_keywords.push_back( "spec.reflect" );
  mcnp5_keywords.push_back( "white.reflect" );
  mcnp5_keywords.push_back( "graveyard" );
  
  mcnp5_keyword_synonyms[ "rest.of.world" ] = "graveyard";
  mcnp5_keyword_synonyms[ "outside.world" ] = "graveyard";

  // parse data from geometry
  rval = DAG->parse_properties( mcnp5_keywords, mcnp5_keyword_synonyms );
  if (MB_SUCCESS != rval) {
    std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  std::string lfname(lfile, *llen);
  std::cerr << "Going to write an lcad file = " << lfname << std::endl;
  // Before opening file for writing, check for an existing file
  if( lfname != "lcad" ){
    // Do not overwrite a lcad file if it already exists, except if it has the default name "lcad"
    if( access( lfname.c_str(), R_OK ) == 0 ){
      std::cout << "DagMC: reading from existing lcad file " << lfname << std::endl;
      return; 
    }
  }

  std::ofstream lcadfile( lfname.c_str() );

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
  }

  // cells finished, skip a line
  lcadfile << std::endl;
  
  // write the surface cards
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

