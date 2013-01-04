#include "moab/Interface.hpp"
#include <iostream>

/* Do ray fire
 * *ih  - Volume ID to do ray fire against
 * *jsu - ? (RefFace ID)
 * *nps - ?
 * (*uuu,*vvv,*www) - Ray direction vector
 * (*xxx,*yyy,*zzz) - Ray point
 * *huge - ?
 * *dls - output distnace of intersection
 * *jap - Next intersected surface, or zero if none
 */
  // void dagmctrack_(int *ih, double *uuu,double *vvv,double *www,double *xxx,
    //                double *yyy,double *zzz,double *huge,double *dls,int *jap,int *jsu,
      //             int *nps );


using namespace moab;


/**
*  fludagW.cpp
*  jcz
*  This is a wrapper, hence the 'W'
*  Implement FLUKA method stubs (originally fortranxx) in c++;  
*  Use print statements to determine how the FLUKA calls work
*/
/*
* Called from 	GEODBG, which executes the geometry debugging.  
*/
extern "C" int LOOKDB (double X, double Y, double Z, double MXERLM)
{
	return 1;
}

/*
* Called from GEOMAG, returns the combinatorial geometry zone of point (X,Y,Z) to determine 
* the region of the end step in magnetic field.
*/
extern "C" int LOOKMG (double X, double Y, double Z,
                  double dir[3], // Direction cosines vector
		  int RegionNum, // region number
                  int *newCell,    // region # of p'le after step ("IRPRIM" in FLUKA)
                  int *Ierr        // error code
)
{
	return 1;
}

/*
 * Tracking routine for solving tracking problems.
 * Same input and output as LOOKMG.
*/
extern "C" int LOOKFX (double X, double Y, double Z,
                  double dir[3], // Direction cosines vector
		  int RegionNum, // region number
                  int *newCell,    // region # of p'le after step ("IRPRIM" in FLUKA)
                  int *Ierr        // error code
)
{
	return 1;
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


// *ih              - volume index
// *uuu, *vvv, *www - ray direction
// *xxx, *yyy, *zzz - ray point
// *huge            - passed to ray_fire as 'huge'
// *dls             - output from ray_fire as 'dist_traveled'
// *jap             - intersected surface index, or zero if none
// *jsu             - previous surface index
/*
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

  // detect streaming or reflecting situations 
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
*/

/*
int main(int argc, char* argv[])
{
	double p0=1;
	double p1=2;
	double p2=3;

	double dir0=5;
	double dir1=5;
	double dir2=0;

	int RegionNum = 5;
	double dist =.01;
	int SurfaceNum = 4;
	double step = 0;
	int newCell = 0;
	double tol  = 0;
	int Ierr = 1;


	double pos[3] = {p0, p1, p2};
	double dir[3] = {dir0, dir1, dir2};

// 	std::cout << " Cartesian coordinate vector: pos[3] = " << pos[0]<<  pos[1]  << std::endl;
	
	int ret = G1FLU_(pos, dir, RegionNum, dist, &SurfaceNum, &step, &newCell, &tol, &Ierr);

	return 1;
}
*/
