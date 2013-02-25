
////////////////////////////////////////////////////////////////////
//
// WrapNorml.hh - Sara Vanini 
//
// Wrapper for computing normal unit-vector in global coordinates.
//
// Fluka requires normal vector exiting from final position (of 
// particle) volume, that is: entering in volume of initial position. 
// Geant4 always computes the normal vector exiting from the volume. 
// In GetLocalExitNormal() call the volume is the pre-step volume 
// (so G4 normal vector sign is opposite of fluka-required normal).
// If IsLocalExitNormalValid=false, normal value is computed from 
// init-step volume (in this case sign must change), or from
// end-step volume (the sign is the same). Normal vector is computed 
// always on boundary between volumes, in global coordinates (to take
// rotation of parameterised volumes in hierarchy in consideration). 
// So: nrmlwr returns inwards pointing unit normal of the shape for 
// surface closest to the point returned by the navigator (last step 
// end-point).
// 
//
////////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"


using namespace moab;
#define DAG DagMC::instance()

void nrmlwr(double& pSx, double& pSy, double& pSz,
            double& pVx, double& pVy, double& pVz,
	    double* norml, const int& oldReg, 
	    const int& newReg, int& flagErr)
{
  std::cout << "============ NRMLWR-DBG =============" << std::endl;
  
  //dummy variables
  flagErr=0;
  
  //return normal:
  norml[0]=0.0;
  norml[1]=0.0;
  norml[2]=0.0;
  std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
}


