#ifndef DAGMC_MCNP_IFACE_H
#define DAGMC_MCNP_IFACE_H

#include <string>
#include <vector>

#include "moab/Types.hpp"

#include "MBInterface.hpp"
#include "MBCartVect.hpp"

#include "DagMC.hpp"

#ifdef __cplusplus
extern "C" {
#endif


/* initialize DAGMC from FORTRAN main 
 * @param max_pbl - The maximum index of the pblcm (temporary particle state) array
 *                  This is the largest n that will arrive in calls to savpar and getpar
 */
  void dagmcinit_(char *cfile, int *clen,  
                  char *ftol,  int *ftlen, 
                  int *parallel_file_mode,
                  double* dagmc_version, int* moab_version, int* max_pbl );

/* 
 * This signature is intended to be called from a c++ function.  
 * It assumes the caller does not need to know the dagmc_version or the moab_version.
 * The original dagmcinit_, which is to be called by a Fortran method, has been reworked 
 * to set up variables and then call this version.
 * 
 * 14 Feb 2013  jcz added boolean running_with_fluka, which prepends "../" to the file name
 */
  void cpp_dagmcinit(const std::string infile, 
                int parallel_file_mode, // parallel read mode
                int max_pbl);

  void slow_check(double pos[3], const double dir[3], int &oldReg);
  MBEntityHandle check_reg(MBEntityHandle volume, double point[3], double dir[3]); // check we are where we say we are
  void special_check(double pos[3],const double dir[3], int& oldReg);

  /* get the sense of a region with respect to the global next_surf,
   * which is set by a call to rayfire
  */
  int getSense(int region);
/*
 * Prepare a descriptive string that creates the properties of the volume whose index is index
 */
  std::string mat_property_string (int index, std::vector<std::string> &properties);
/*
 * Write the material assignment for each volume to a file named matfile
 */
  void fludagwrite_assignma(std::string matfile);  
  void fludagwrite_mat(std::string fname);

/* Add the current particle state to the bank */
  void dagmc_bank_push_( int* nbnk );

/* Revert to the most recently banked particle state */
  void dagmc_bank_usetop_( ) ;
  
/* Remove and forget the most recently banked particle state.
 * Called after dagmc_bank_usetop_() */
  void dagmc_bank_pop_( int* nbnk );

/* Remove all entries in the particle bank */
  void dagmc_bank_clear_( ) ; 

/* Save the current particle state to temporary index n */
  void dagmc_savpar_( int* n );

/* Reset the current particle state using temporary index n*/
  void dagmc_getpar_( int* n );

/* write facet file after initialization and OBBTree generation */
  void dagmcwritefacets_(char *ffile, int *flen);

/* parse metadata and write applications specific data for: MCNP5 */
  void dagmcwritemcnp_(char *lfile, int *llen);

/* Get normal of surface with id *jsu at location (*xxx,*yyy,*zzz) and store
   in three doubles at ang (an arry of length 3) */
  void dagmcangl_(int *jsu, double *xxx, double *yyy, double *zzz, double *ang);


  /* Given point and direction, determine if the particle is going into 
   * cell i1.  Assume the point is on the surface jsu.  Return j=0 if
   * the particle is directed into the cell, j=1 if it is leaving the cell.
   * This is used as a point-in-volume query for points that are known
   * to be on a surface.
   */
  void dagmcchkcel_by_angle_( double *uuu, double *vvv, double *www, 
                              double *xxx, double *yyy, double *zzz,
                              int *jsu, int *i1, int *j);

  /* Point-in-volume query.  Determine if the particle at given coordinates
   * is inside or outside of cell i1.  Return j=1 if outside or on boundary,
   * and j=0 if inside.  
   */
  void dagmcchkcel_(double *uuu,double *vvv,double *www,double *xxx,
                    double *yyy,double *zzz, int *i1, int *j);

/* Determine distance to nearest surface
 * *ih - current RefVolume ID
 * *xxx, *yyy, *zzz - current point
 * *huge - passed definition of a large number
 * *dbmin - Output, distance to nearest surface
 */
  void dagmcdbmin_( int *ih, 
                  double *xxx, double *yyy, double *zzz, 
                  double *huge, double* dbmin);


/* Get next Cell
 * *jsu - Surface to cross
 * *icl - Previous cell ID
 * *iap - Next cell ID (output parameter)
 */
  void dagmcnewcel_( int *jsu, int *icl, int* iap );

  /**
   * Tell dagmc that a particle has changed direction at the most recently reached surface boundary,
   * and reset ray history as necessary.
   * The parameters are the new particle direction and a flag.
   * If *verify_dir_change==0, assume without checking that a direction change has actually happened.
   * If *verify_dir_change==1, compare the new uvw to the last known one, and hold on to 
   * all ray history if not.  (This mode is used during electron transport,
   * which occasionally calls this function without having changed the direction;
   * it's less invasive to do the check in C++ instead of Fortran.)
   */
  void dagmc_surf_reflection_( double *uuu, double *vvv, double *www, int* verify_dir_change );

  void dagmc_particle_terminate_( );

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
  void dagmctrack_(int *ih, double *uuu,double *vvv,double *www,double *xxx,
                   double *yyy,double *zzz,double *huge,double *dls,int *jap,int *jsu,
                   int *nps );

/* Measure entities
 * vols - 2xN array where first column contains, as output, measure of every volume.
 * aras - 2xN array where first column contains, as output, measure of every surface
 */                        
  void dagmcvolume_(int* mxa, double* vols, int* mxj, double* aras);

/* Set distance limit */
  void dagmc_setdis_(double *d);

  void dagmc_set_settings_(int* use_dist_limit, int* use_cad, double* overlap_thickness, int* srccell_mode);

  void dagmc_init_settings_(int* use_dist_limit, int* use_cad,     
                            double* overlap_thickness, double* facet_tol, int* srccell_mode );

  void dagmc_version_(double* dagmcVersion);
#ifdef __cplusplus
} // extern "C"
#endif

#endif /* DAGMC_MCNP_IFACE_H */
