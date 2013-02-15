// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapLookZ.hh - Sara Vanini 
//
// Wrapper for localisation of starting point of particle.
//
// modified 20/III/00: history initialization moved to ISVHWR
// modified 13/IV/00: located history saved in jrLtcGeant and 
//                    incremented
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
// modified 17/06/02: by I. Gonzalez. STL migration
// modified 27/11/09: by P.Sala added compatibility with Plotgeom
//
//////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
// #include "FGeometryInit.hh"
/*
#include "G4VPhysicalVolume.hh"
#include "FluggNavigator.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalVolumeStore.hh"
*/
// #include "globals.hh"


void lkwr(double& pSx, double& pSy, double& pSz,
          double* pV, const int& oldReg, const int& oldLttc,
          int& newReg, int& flagErr, int& newLttc)
{
  //flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "======= LKWR =======" << std::endl;
#endif
 /* 
  //FGeometryInit, navigator, volumeStore  pointers
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  FluggNavigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();

  int geoflg =   ptrGeoInit-> getGeoFlag();
  //coordinates in mm.
  G4ThreeVector pSource(pSx,pSy,pSz);
  pSource *= 10.0; //in millimeters!
  G4ThreeVector pVec(pV[0],pV[1],pV[2]);
    
  //locate point and update histories
  G4TouchableHistory * ptrTouchableHistory = 
    ptrGeoInit->GetTouchableHistory();
  ptrNavig->LocateGlobalPointAndUpdateTouchable(pSource,pVec,
						ptrTouchableHistory,true);
  //updating tmp but not old histories, they are useful in 
  //case of RGRPWR call, or when fluka, after a LOOKZ call, 
  //descards step for multiple scattering and returns to old history
  //NO, after lattice-fix we don't need old history anymore!
  
  ptrGeoInit->UpdateHistories(ptrTouchableHistory->GetHistory(),0);
  G4VPhysicalVolume * located = ptrTouchableHistory->GetVolume();
  
  //if volume not found, out of mother volume: returns "number of volumes"+1
  if(!located) {
    int numVol = int(pVolStore->size());
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "Out of mother volume! " << numVol << std::endl ;
#endif
    newReg = numVol + 1;
  }
  else { 
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* ISVHWR call to store current NavHistWithCount in jrLtGeant"
	   << std::endl;
#endif 
    
    //save history in jrLtGeant and increment counter  
    int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
    int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
    // dummy flag for plotgeom
    if (geoflg != 2003 )
      {
     LttcFlagGeant += 1;
    jrLtGeant[LttcFlagGeant] = isvhwr(0,0);
    
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* CONHWR call to increment counter" << std::endl;
#endif 
    int incrCount=1;
    conhwr(jrLtGeant[LttcFlagGeant],&incrCount);
    
    //update LttcFlagGeant
    ptrGeoInit->SetLttcFlagGeant(LttcFlagGeant);
      }
    //return region number and dummy variables
    int volIndex = ~0;
    for (unsigned int i=0; i<pVolStore->size(); i++)
      if ((*pVolStore)[i] == located) volIndex = i;
    //int volIndex=int(pVolStore->index(located));
    if (volIndex==(~0)) {
      std::cerr << "FLUGG: Problem in routine WrapG1 tryingto find volume after step" << std::endl;
      exit(-999);
    }
    newReg=volIndex+1;   
    if (geoflg != 2003 )    newLttc=jrLtGeant[LttcFlagGeant];
    flagErr=newReg;
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "LKWR Located Physical volume = ";
    std::cout << located->GetName() << std::endl;
#endif
  }
    flagErr=newReg;
*/
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "======= Out of LKWR =======" << std::endl;
#endif
}






