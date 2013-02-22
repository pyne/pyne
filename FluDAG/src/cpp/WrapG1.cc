

///////////////////////////////////////////////////////////////////
//
// WrapG1.hh - Sara Vanini
//
// Wrapper for geometry tracking: returns approved step of 
// particle and alla variables that fluka G1 computes.
//
//
/////////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

void g1wr(double& pSx, double& pSy, double& pSz, double* pV,
          int& oldReg, const int& oldLttc, double& propStep,
          int& nascFlag, double& retStep, int& newReg,
          double& saf, int& newLttc, int& LttcFlag,
          double* sLt, int* jrLt)
{
  std::cout<<"============= G1WR =============="<<std::endl;    
}




