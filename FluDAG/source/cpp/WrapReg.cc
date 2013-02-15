
// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapReg.hh - Sara Vanini 
//
// Wrapper for scoring hits: previous step end-point is taken from 
// history (and compared with fluka region index, flukaReg),
// then the wrapper returns all the information regarding the 
// volume tree, i.e. returns indMother[] array with all the 
// mother volumes index and repMother[] array with all the 
// mother volumes repetition number.   
//
// modified: 16/III/99
// modified: 14/IV/00 ptrLttc included 
// modified: 24.10.00: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
// modified: 17/06/02 by I. Gonzalez. STL migration.
//
///////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
// #include "FGeometryInit.hh"
/*
#include "NavHistWithCount.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalVolumeStore.hh"
*/
// #include "globals.hh"


void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
            int* indMother, int* repMother, int& depthFluka)
{
  //flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "============= RGRPWR ==============" << std::endl;    
  std::cout << "ptrLttc=" << ptrLttc << std::endl;
#endif 
  
/*
  //Geoinit, Navigator, VolStore pointers
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  
  //get jrLtGeant array and flag
  int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
  int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
  
  G4bool foundHistory = false;
  int i = LttcFlagGeant;
  
  while(!foundHistory && i>=0) {
    if(jrLtGeant[i]==ptrLttc) foundHistory = true;
    i -= 1;
  }
  
  if(!foundHistory) {
    std::cout << "* ERROR! History not in jrLtGeant!" << std::endl;
    //only in debugging version....
    assert(foundHistory);
  }
  else {
    //get history pointer from ptrLttc
    NavHistWithCount* ptrNavHistCount = 
      reinterpret_cast<NavHistWithCount*>(ptrLttc);
    G4NavigationHistory* ptrNavHist = ptrNavHistCount->GetNavHistPtr();
    
    G4VPhysicalVolume* ptrVolHistory = 0;
    int volHistIndex = 0;
    int depth = ptrNavHist->GetDepth();
    for(int h=0; h<=depth; h++) {
      ptrVolHistory = ptrNavHist->GetVolume(h);
      //
      volHistIndex = ~0;
      for (unsigned int i=0; i<pVolStore->size(); i++)
	if ((*pVolStore)[i] == ptrVolHistory) 
	  volHistIndex = i; 
      if (volHistIndex==(~0)) {
	G4cerr << "FLUGG: Problem in routine WrapReg tryingto find volume after step" << std::endl;
	exit(-999);
      }
      //volHistIndex = int(pVolStore->index(ptrVolHistory));
      //
      indMother[h] = volHistIndex+1;
      if(ptrVolHistory->IsReplicated()) {
	//true if volume is replica or parameterized,
	//false for placements; repetition numbers 
	//are set: 1,2,3,etc.; 0 for placed volumes.
	repMother[h] = 1+int(ptrNavHist->GetReplicaNo(h));
      }
      else {
	repMother[h] = 0;
      }
      
#ifdef  DAGGEOMETRY_DEBUG
      //  std::cout<<"Level="<<h<<"   : ";
      //  std::cout<<"Region="<<indMother[h]<<"   (repetition="<<
      //    repMother[h]<<")"<<std::endl;
#endif
    }
    
    //compute new region index
    int volIndex = ~0;
    for (unsigned int i=0; i<pVolStore->size(); i++)
      if ((*pVolStore)[i] == ptrVolHistory) 
	volIndex = i;
    //int volIndex=int(pVolStore->index(ptrVolHistory));
    if (volIndex==(~0)) {
      G4cerr << "FLUGG: Problem in routine WrapReg tryingto find volume after step" << std::endl;
      exit(-999);
    }
    
    g4Reg=volIndex+1;
    
    depthFluka = depth;
  }
 */ 
}





