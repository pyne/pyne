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

  void dagmc_version_(double* dagmcVersion);
#ifdef __cplusplus
} // extern "C"
#endif

#endif /* DAGMC_MCNP_IFACE_H */
