


///////////////////////////////////////////////////////////////////
//
// Interface for Dag Wrappers
//
///////////////////////////////////////////////////////////////////

#ifndef DAGWRAPPERS_HH
#define DAGWRAPPERS_HH

// #include "DagMC.hpp"
// #define DAG DagMC::instance()
// #include "globals.hh"

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

// Wrapflgfwr

extern "C" void flgfwr(int& flkflg);
// WrapDN

extern "C" int idnrwr(const int & nreg, const int & mlat);

// WrapG1

extern "C" void  g1wr(double& pSx, double& pSy, double& pSz, double* pV,
                      int& oldReg, const int& oldLttc, double& propStep,
                      int& nascFlag, double& retStep, int& newReg,
	              double& saf, int& newLttc, int& LttcFlag,
                      double* sLt, int* jrLt);

// WrapG1RT

extern "C" void g1rtwr(void);

// WrapIncrHist

extern "C" void conhwr(int& intHist, int* incrCount); 

// WrapIniHist

extern "C" void inihwr(int& intHist);                   

// WrapInit - calls createFlukaMatRegion(...)

extern "C" void jomiwr(int & nge, const int& lin, const int& lou,
                       int& flukaReg);

// WrapLookDB

extern "C" void lkdbwr(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
               	       int& newReg, int& flagErr, int& newLttc);

// WrapLookFX

extern "C" void lkfxwr(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
                       int& newReg, int& flagErr, int& newLttc);
	    
// WrapLookMG

extern "C" void lkmgwr(double& pSx, double& pSy, double& pSz,
                       double* pV, const int& oldReg, const int& oldLttc,
		       int& flagErr, int& newReg, int& newLttc);
	    
// WrapLookZ

extern "C" void lkwr(double& pSx, double& pSy, double& pSz,
                     double* pV, const int& oldReg, const int& oldLttc,
	             int& newReg, int& flagErr, int& newLttc);

// WrapMag

extern "C" void fldwr(const double& pX, const double& pY, const double& pZ,
                       double& cosBx, double& cosBy, double& cosBz, 
                       double& Bmag, int& reg, int& idiscflag);
	    
// WrapNorml

extern "C" void nrmlwr(double& pSx, double& pSy, double& pSz,
                       double& pVx, double& pVy, double& pVz,
	               double* norml, const int& oldReg, 
	               const int& newReg, int& flagErr);

// WrapReg

extern "C" void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
                       int* indMother, int* repMother, int& depthFluka);

// WrapSavHist
	    
extern "C" int isvhwr(const int& fCheck, const int& intHist);

//Wrap region2name

extern "C" void rg2nwr(const int& mreg, char* Vname);

#endif //WRAPPERS_HH

