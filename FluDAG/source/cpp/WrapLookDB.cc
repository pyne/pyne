
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapG1.hh - Sara Vanini
//
// Dummy wrapper (in fluka: for geometry debugging)
//
//////////////////////////////////////////////////////////////////

//#ifndef lkdbwr
//#define lkdbwr lkdbwr_

#include "DagWrappers.hh"
// #include "globals.hh"

void lkdbwr(double& pSx, double& pSy, double& pSz,
	    double* pV, const int& oldReg, const int& oldLttc,
	    int& newReg, int& flagErr, int& newLttc)
{
  //flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout<<"============= LKDBWR =============="<< std::endl;
#endif
  
  //return region number and dummy variables
  newReg=0;   
  newLttc=0;
  flagErr=-1; 
}
//#endif
