
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
// modified 10/III/99
// modified 25/V/00
// modified 7/VI/00 for boundary-crossing in case of relocation
// modified 5/VII/00 geometry error on boundary fixed
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
////////////////////////////////////////////////////////////////////

/*
#include "Wrappers.hh"
#include "WrapUtils.hh"
#include "FGeometryInit.hh"
#include "G4VPhysicalVolume.hh"
#include "FluggNavigator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
*/
#include "DagWrappers.hh"
#include "DagWrapUtils.hh"


using namespace moab;
#define DAG DagMC::instance()

void nrmlwr(double& pSx, double& pSy, double& pSz,
            double& pVx, double& pVy, double& pVz,
	    double* norml, const int& oldReg, 
	    const int& newReg, int& flagErr)
{
  //flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "============ NRMLWR-DBG =============" << std::endl;
#endif
  
  //dummy variables
  flagErr=0;
  
  //navigator pointer
  // static FGeometryInit * ptrGeoInit;
  // ptrGeoInit = FGeometryInit::GetInstance();
  // FluggNavigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
 /* 
  //variables
  G4ThreeVector normalLoc;
  G4ThreeVector normalGlob;
  
  //normal computing
  //  if(ptrNavig->IsExitNormalValid()) {
  //  normalLoc=ptrNavig->GetLocalExitNormal();
  //  normalLoc *= -1;        
    
  //  //global cooordinates normal
  //  normalGlob = ptrNavig->GetLocalToGlobalTransform().
  //    TransformAxis(normalLoc);
  // }
  G4bool valid;
  normalLoc=ptrNavig->GetLocalExitNormal(&valid);
  if(valid) {
    normalLoc *= -1;        
    
    //global cooordinates normal
    normalGlob = ptrNavig->GetLocalToGlobalTransform().
      TransformAxis(normalLoc);
  }
  else {
    G4VPhysicalVolume *touchvolume;
    G4ThreeVector theLocalPoint;
    
    if(ptrNavig->EnteredDaughterVolume()) {
      //volume from touchable history
      G4TouchableHistory *ptrTouchableHistory=ptrGeoInit->
	GetTouchableHistory();
      touchvolume = ptrTouchableHistory->GetVolume();
      
      //local point from navigator and normal
      theLocalPoint = ptrNavig->GetCurrentLocalCoordinate(); 
      normalLoc = touchvolume->GetLogicalVolume()->GetSolid()->
	SurfaceNormal(theLocalPoint);
      
      //global cooordinates normal
      normalGlob = ptrNavig->GetLocalToGlobalTransform().
	TransformAxis(normalLoc);
    }
    else {
      //volume from old history
      const G4NavigationHistory * ptrOldNavHist = 
	ptrGeoInit->GetOldNavHist()->GetHistory();
      touchvolume = ptrOldNavHist->GetTopVolume();
      
      if(!touchvolume) {
	// Old history has been reseted by LOOKZ relocation,
	// so is necessary to track back and forward to find
	// the right histories.
	
	
	////////////  COMPUTE STEP BACKWARD //////////////////
	G4ThreeVector theGlobalPoint = ptrNavig->GetLocalToGlobalTransform().
	  TransformPoint(ptrNavig->GetCurrentLocalCoordinate());
	
	//compute step and new location
	G4ThreeVector pVec(pVx,pVy,pVz);
	pVec=-pVec;
	double appStep = 0;
	double safety = 0;
	G4bool onBoundary = false;
	double physStep = 1000000000;
	int newRegStep;
	G4ThreeVector partLoc =  theGlobalPoint;
	G4bool fErr=false;
	
#ifdef DAGGEOMETRY_DEBUG
	std::cout << "Old history not found" << std::endl;
	std::cout << "* NRML needs boundary-crossing: computing step backward..."
	       << std::endl;
#endif
	
	//compute step and location 
	newRegStep=StepAndLocation(partLoc,pVec,physStep,
				   appStep,safety,onBoundary,
				   fErr,oldReg);
	
	G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::
	  GetInstance();
	
	if(appStep<physStep && newRegStep!=int(pVolStore->size())+1) {
	  //end-step point is on boundary between volumes;
	  //==> update step-histories
	  
#ifdef DAGGEOMETRY_DEBUG
	  std::cout << "* updating step-histories" << std::endl;
#endif
	  
	  G4TouchableHistory * ptrTouchHist = 
	    ptrGeoInit->GetTouchableHistory();
	  ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),2);
	}
	else {
#ifdef DAGGEOMETRY_DEBUG
	  std::cout << "ERROR! Boundary not found" << std::endl;
#endif
	}
	
#ifdef DAGGEOMETRY_DEBUG
	std::cout << "* computing step forward..." << std::endl;
#endif
	pVec=-pVec;
	safety = 0;
	onBoundary = false;
	
	//compute step and location for boundary crossing
	newRegStep=StepAndLocation(partLoc,pVec,physStep,
				   appStep,safety,onBoundary,fErr,oldReg);
	if(appStep<physStep) {
	  //end-step point is on boundary between volumes;
	  //==> update step-histories
	  
#ifdef DAGGEOMETRY_DEBUG
	  std::cout << "* updating step-histories" << std::endl;
#endif
	  
	  G4TouchableHistory * ptrTouchHist = 
	    ptrGeoInit->GetTouchableHistory();
	  ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),2);
	}
      }
      
      // now touchvolume exist.
      // local point from navigator and global point
      // N.B. if particle has exited world volume, 
      // FluggNavigator doesn't update lastLocatedPoint. 
      // So be carefull in building geometry always to have a 
      // big world volume that fluka won't exit.
      
      touchvolume = ptrOldNavHist->GetTopVolume();
      
      G4ThreeVector theGlobalPoint = ptrNavig->GetLocalToGlobalTransform().
	TransformPoint(ptrNavig->GetCurrentLocalCoordinate());
      theLocalPoint = ptrOldNavHist->GetTopTransform().
	TransformPoint(theGlobalPoint);
      normalLoc = (touchvolume->GetLogicalVolume()->GetSolid()->
		   SurfaceNormal(theLocalPoint));
      normalLoc *= -1; 
      
      //global cooordinates normal
      normalGlob = ptrOldNavHist->GetTopTransform().
	Inverse().TransformAxis(normalLoc);	    
    }
  }
  
  //return normal:
  norml[0]=double(normalGlob.x());
  norml[1]=double(normalGlob.y());
  norml[2]=double(normalGlob.z());
*/  
  //return normal:
  norml[0]=0.0;
  norml[1]=0.0;
  norml[2]=0.0;
  //for debugging
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "Normal: " << norml[0] << ", " << norml[1] << ", " << norml[2]  << std::endl;
#endif
}


