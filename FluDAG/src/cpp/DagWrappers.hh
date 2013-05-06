

///////////////////////////////////////////////////////////////////
//
// Interface for Dag Wrappers
//
///////////////////////////////////////////////////////////////////

#ifndef DAGWRAPPERS_HH
#define DAGWRAPPERS_HH


#define idnrwr idnrwr_
#define g1wr g1wr_
#define g1rtwr g1rtwr_
#define conhwr conhwr_
#define inihwr inihwr_
#define jomiwr jomiwr_
#define lkdbwr lkdbwr_
#define lkfxwr lkfxwr_
#define lkmgwr lkmgwr_
#define lkwr lkwr_
#define fldwr fldwr_
#define flgfwr flgfwr_
#define nrmlwr nrmlwr_
#define rgrpwr rgrpwr_
#define isvhwr isvhwr_
#define rg2nwr rg2nwr_

// WrapFlgfwr.cc stubs this
extern "C" void flgfwr(int& flkflg);

// WrapDN.cc has been deleted.  The FLUKA internal function is used.
extern "C" int idnrwr(const int & nreg, const int & mlat);

// WrapG1.cpp has been deleted.  The function is now
// defined in fluka_funcs.cpp.  It calls g1_fire.
extern "C" void  g1wr(double& pSx, double& pSy, double& pSz, double* pV,
                      int& oldReg, const int& oldLttc, double& propStep,
                      int& nascFlag, double& retStep, int& newReg,
	              double& saf, int& newLttc, int& LttcFlag,
                      double* sLt, int* jrLt);

  /**
   * g1_fire is called by fludag's implementation of g1wr.  It calls DAG->ray_fire(...).
   * oldRegion region of start point
   * point     start point
   * dir       direction vector
   * propStep
   * retStep   set to distance to next surface
   * newRegion region ofter step
  */
 //  void g1_fire(int& oldRegion, double point[], double dir[],  double& retStep, int& newRegion);
  void g1_fire(int& oldRegion, double point[], double dir[], 
               double &propStep, double& retStep,  int& newRegion);


// WrapG1RT.cc has been deleted.  The function is now 
// defined in fluka_funcs.cpp, as a stub.
extern "C" void g1rtwr(void);

// WrapIncrHist.cc has been deleted.  The internal 
// FLUKA form of this function is called.
// extern "C" void conhwr(int& intHist, int* incrCount); 

// WrapIniHist.cc has been deleted.  The function is now
// commented out in fluka_funcs.cpp
// extern "C" void inihwr(int& intHist);                   

// WrapInit.cc - has been deleted, the function is now
// defined in fluka_funcs.cpp
extern "C" void jomiwr(int & nge, const int& lin, const int& lou,
                       int& flukaReg);

// WrapLookDB.cc has been deleted, the function is now
// defined in fluka_funcs.cpp.  It sets some of its return values
extern "C" void lkdbwr(double& pSx, double& pSy, double& pSz,
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
extern "C" void lkwr(double& pSx, double& pSy, double& pSz,
                     double* pV, const int& oldReg, const int& oldLttc,
	             int& newReg, int& flagErr, int& newLttc);

// Wrapper for lkwr clarifying which areguments are used.
int look( double& posx, double& posy, double& posz, double* dir, int& region);

// WrapMag.cc stubs this function and is linked in
extern "C" void fldwr(const double& pX, const double& pY, const double& pZ,
                       double& cosBx, double& cosBy, double& cosBz, 
                       double& Bmag, int& reg, int& idiscflag);
	    
// WrapNorml.cc has been removed.  This function is defined in fluka_funcs.
//  It is called by normal(..).
extern "C" void nrmlwr(double& pSx, double& pSy, double& pSz,
                       double& pVx, double& pVy, double& pVz,
	               double* norml, const int& oldReg, 
	               const int& newReg, int& flagErr);

// Wrapper for nrmlwr clarifying which arguments are used
int  normal (double& posx, double& posy, double& posz, double *norml, int& regionForSense);

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

#endif //DAGWRAPPERS_HH

