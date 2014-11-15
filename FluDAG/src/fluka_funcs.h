#ifndef FLUDAG_SRC_FLUKA_FUNCS_H
#define FLUDAG_SRC_FLUKA_FUNCS_H

#include <iostream>
#include <stdlib.h>
#include <string>      // string, cout
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
  // check we are where we say we are
  MBEntityHandle check_reg(MBEntityHandle volume, double point[3], double dir[3]); 

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


///////////////////////////////////////////////////////////////////
//
// Interface for Dag Wrappers
//
///////////////////////////////////////////////////////////////////

// #include "MBInterface.hpp"
// #include "MBCartVect.hpp"
// #include "moab/Types.hpp"

#define f_idnr idnrwr_
#define g_step g1wr_
#define f_g1rt g1rtwr_
#define inihwr inihwr_
#define jomiwr jomiwr_
#define f_lookdb lkdbwr_
#define lkfxwr lkfxwr_
#define lkmgwr lkmgwr_
#define f_look lkwr_
#define fldwr fldwr_
#define flgfwr flgfwr_
#define f_normal nrmlwr_
#define rgrpwr rgrpwr_
#define isvhwr isvhwr_
#define rg2nwr rg2nwr_

// WrapFlgfwr.cc stubs this
extern "C" void flgfwr(int& flkflg);

// The FLUKA internal function is used.
extern "C" int f_idnr(const int & nreg, const int & mlat);

// The function is defined in fluka_funcs.cpp.  It calls g_fire.
extern "C" void  g_step(double& pSx, double& pSy, double& pSz, double* pV,
                      int& oldReg, const int& oldLttc, double& propStep,
                      int& nascFlag, double& retStep, int& newReg,
	              double& saf, int& newLttc, int& LttcFlag,
                      double* sLt, int* jrLt);

/**
 * g_fire is called by fludag's implementation of g_step.  It calls DAG->ray_fire(...).
 * oldRegion region of start point
 * point     start point
 * dir       direction vector
 * propStep
 * retStep   set to distance to next surface
 * newRegion region ofter step
 */
void g_fire(int& oldRegion, double point[], double dir[], 
	      double &propStep, double& retStep, double &safety,
	      int& newRegion);

// Stub function
extern "C" void f_g1rt(void);

// WrapInit.cc - has been deleted, the function is now
// defined in fluka_funcs.cpp
extern "C" void jomiwr(int & nge, const int& lin, const int& lou,
                       int& flukaReg);

// WrapLookDB.cc has been deleted, the function is now
// defined in fluka_funcs.cpp.  It sets some of its return values
extern "C" void f_lookdb(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
               	       int& newReg, int& flagErr, int& newLttc);

// WrapLookFX
// Stubbed in WrapLookFX.cc and linked in.
extern "C" void lkfxwr(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
                       int& newReg, int& flagErr, int& newLttc);
	    
// WrapLookMG.cc stubs this function and is linked in.
extern "C" void lkmgwr(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
		       int& flagErr, int& newReg, int& newLttc);
	    
// WrapLookZ has been deleted.  This function is defined in
// fluka_funcs.cc.  It is called by look(..).
extern "C" void f_look(double& pSx, double& pSy, double& pSz,
                     double* pV, const int& oldReg, const int& oldLttc,
	             int& newReg, int& flagErr, int& newLttc);

// Wrapper for f_look clarifying which arguments are used.
int look( double& posx, double& posy, double& posz, double* dir, int& region);

// WrapMag.cc stubs this function and is linked in
extern "C" void fldwr(const double& pX, const double& pY, const double& pZ,
                       double& cosBx, double& cosBy, double& cosBz, 
                       double& Bmag, int& reg, int& idiscflag);
	    
// WrapNorml.cc has been removed.  This function is defined in fluka_funcs.
//  It is called by normal(..).
extern "C" void f_normal(double& pSx, double& pSy, double& pSz,
                       double& pVx, double& pVy, double& pVz,
	               double* norml, const int& oldReg, 
	               const int& newReg, int& flagErr);

// Wrapper for f_normal clarifying which arguments are used
int  normal (double& posx, double& posy, double& posz, double *norml, int& regionForSense);

// Protoype for boundary test funct
int boundary_test(MBEntityHandle vol, double xyz[3], double uvw[3]);

// protoype for mat props
std::string mat_property_string (int index, std::vector<std::string> &properties);

// WrapReg

extern "C" void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
                       int* indMother, int* repMother, int& depthFluka);

// WrapSavHist.cc has been removed.  This function has been commented
// out of an implementation in fluka_funcs.cpp, so the internal FLUKA 
// version, if any, is called.
extern "C" int isvhwr(const int& fCheck, const int& intHist);

// WrapReg2name.cc defines this and is linked in.  It calls
// region2name, which is defined in fluka_funcs
// ToDo:  Remove WrapReg2name.cc and implement rg2nwr(..) in fluka_funcs.cpp
extern "C" void rg2nwr(const int& mreg, char* Vname);
// Defined in fluka_funcs, called by rg2nwr
void region2name(int volindex, char * vname );


#endif /* FLUDAG_SRC_FLUKA_FUNCS_H */
