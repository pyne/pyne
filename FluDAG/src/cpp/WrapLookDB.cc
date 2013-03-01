
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapG1.hh - Sara Vanini
//
// Dummy wrapper (in fluka: for geometry debugging)
//
//////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

void lkdbwr(double& pSx, double& pSy, double& pSz,
	    double* pV, const int& oldReg, const int& oldLttc,
	    int& newReg, int& flagErr, int& newLttc)
{
  std::cerr<<"============= LKDBWR =============="<< std::endl;
  
  //return region number and dummy variables
  newReg=0;   
  newLttc=0;
  flagErr=-1; 

  return;
}
