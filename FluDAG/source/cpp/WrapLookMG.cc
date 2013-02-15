
// Dag tag 

///////////////////////////////////////////////////////////////////
//
// WrapLookMG.hh - Sara Vanini 26/X/99
//
// Wrapper for localisation of particle for magnetic field tracking 
//
// modified 13/IV/00: check if point belongs to one of the lattice 
// histories stored in jrLtGeant 
//
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
// modified 17/06/02: by I. Gonzalez. STL migration
//
///////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
// #include "FGeometryInit.hh"
/*
#include "NavHistWithCount.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4NormalNavigation.hh"
#include "G4VoxelNavigation.hh"
#include "G4ParameterisedNavigation.hh"
#include "G4ReplicaNavigation.hh"
*/
// #include "globals.hh"


//auxiliary function declarations
/*
bool PointLocate(const G4NavigationHistory &,
			  const G4ThreeVector &, 
			  const G4ThreeVector *,
			  int &);
EVolume CharacteriseDaughters(const G4LogicalVolume *);
*/

void lkmgwr(double& pSx, double& pSy, double& pSz,
            double* pV, const int& oldReg, const int& oldLttc,
	    int& flagErr, int& newReg, int& newLttc)
{
//flag
#ifdef DAGGEOMETRY_DEBUG
    std::cout<<"============= LKMGWR =============="<<std::endl;
#endif

/*
    //Geoinit, Navigator, etc. pointers
    static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
    G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
    
    //setting variables (and dimension: Fluka uses cm.!)
    G4ThreeVector globalPointcm(pSx,pSy,pSz);
    G4ThreeVector globalPoint =  globalPointcm * 10.; //in mm.
    G4ThreeVector globalDirection(pV[0],pV[1],pV[2]);
    
    //get jrLtGeant array and initialize variables
    int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
    int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
    bool belongToVolume = false;
    int i = LttcFlagGeant;
 //P.S. 6 jul 2009
    int regionIn = oldReg;
    newReg = -1;
    newLttc = -1;
    
    while(!belongToVolume && i>=0) {
      //get history from jrLtGeant
      NavHistWithCount * ptrNavHistCount = 
	reinterpret_cast<NavHistWithCount*>(jrLtGeant[i]);
      const G4NavigationHistory * ptrNavHist =
	ptrNavHistCount->GetNavHistPtr();
      
      //check if globalPoint belongs to volume (can't call 
      //LocateGlobalPoint... because of flag settings)
      belongToVolume = PointLocate(*ptrNavHist,globalPoint,
				   &globalDirection,regionIn);
      
      
      //if point belongs to surface, check scalar product direction*normal
      if(regionIn==-100) {
#ifdef DAGGEOMETRY_DEBUG
	std::cout<<"On surface!"<<std::endl;
#endif
	
	//if entering, then point belongs to volume, 
	//if exiting, point doesn't belong to volume.
	int oldReg,flagErr,reg;
	double x,y,z,px,py,pz;
	double * norml = new double[3];
	x = globalPoint.x();
	y = globalPoint.y();
	z = globalPoint.z();
	px = globalDirection.x();
	py = globalDirection.y();
	pz = globalDirection.z();
	
	nrmlwr(x,y,z,px,py,pz,norml,oldReg,reg,flagErr);
	
	G4ThreeVector normal(norml[0],norml[1],norml[2]);
#ifdef DAGGEOMETRY_DEBUG
	// std::cout<<"Scalar product="<<globalDirection.dot(normal)<<std::endl;
#endif
	if(globalDirection.dot(normal)>0) {
#ifdef DAGGEOMETRY_DEBUG
	  std::cout<<"entering volume!"<<std::endl;
#endif
	  int volIndex = ~0;
	  for (unsigned int i=0; i<pVolStore->size(); i++)
	    if ((*pVolStore)[i] == ptrNavHist->GetTopVolume()) 
	      volIndex = i;
	  if (volIndex==(~0)) {
	    std::cerr << "FLUGG: Problem in routine WrapG1 tryingto find volume after step" << std::endl;
	    exit(-999);
	  }
	  //regionIn = int(pVolStore->index(ptrNavHist->GetTopVolume()))+1;
	  regionIn = volIndex+1; 
	  belongToVolume = true;
	}
      }
      
      i -= 1;
    }
    
    //output variables
    if(belongToVolume) {
      newReg = regionIn;
      newLttc = jrLtGeant[i+1];
    }

    flagErr = pVolStore->size() + 2;      //flagErr=fluka region number + 1
    
    
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "Global point (cm): " << globalPointcm << std::endl;
    std::cout << "Direction: " << globalDirection << std::endl;
    if(newReg!=-1) 
      std::cout << "Point belongs to region " << newReg
	     << " ptr history=" << newLttc << std::endl;
    else 
      std::cout << "No containing volume found!" << std::endl;
#endif
*/
}




//returns false if the point doesn't belong to history top volume (can either 
//belong to one of its daughter volumes or none), otherwise returns true.

/*
bool PointLocate(const G4NavigationHistory & blockedHistConst,
		   const G4ThreeVector & globalPoint, 
		   const G4ThreeVector * pGlobalDirection,
		   int & reg)
{
  
  //variables
  G4NavigationHistory * fHistory = new G4NavigationHistory(blockedHistConst);
  bool belongVolume;
  
  //G4 flags resetted (see: ResetStackAndState)
  bool fExiting = false; 
  G4VPhysicalVolume * fBlockedPhysicalVolume=0;
  int fBlockedReplicaNo=-1;
  bool fLocatedOnEdge=false;  
  
  bool notKnownContained=true,noResult;
  G4VPhysicalVolume *targetPhysical;
  G4LogicalVolume *targetLogical;
  G4VSolid *targetSolid;
  G4ThreeVector localPoint,localDirection;
  EInside insideCode;
  
  
  //Helpers/Utility classes
  G4NormalNavigation  fnormalNav;
  G4VoxelNavigation fvoxelNav;
  G4ParameterisedNavigation fparamNav;
  G4ReplicaNavigation freplicaNav;	  
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  
  //local variables
  localPoint=fHistory->GetTopTransform().TransformPoint(globalPoint);
  localDirection=fHistory->GetTopTransform().TransformPoint(*pGlobalDirection);
  
     //P.S. 6 jul:
    // do not relocate in black hole
    int numVol = int(pVolStore->size());
    if ( reg == numVol+1 ) return true ;
    else
      reg=0;
    //End PS
 //search if volume contains point
  if (fHistory->GetTopVolumeType()!=kReplica) {
    targetSolid=fHistory->GetTopVolume()->GetLogicalVolume()->GetSolid();
    insideCode=targetSolid->Inside(localPoint);
  }
  else {
    insideCode=freplicaNav.BackLocate(*fHistory,globalPoint,
				      localPoint,fExiting,notKnownContained);
    // !CARE! if notKnownContained returns false then the point is within
    // the containing placement volume of the replica(s). If insidecode
    // will result in the history being backed up one level, then the
    // local point returned is the point in the system of this new level
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "Replica: fExiting=" << fExiting << std::endl;
    std::cout << "         notKnownContained=" << notKnownContained << std::endl;
#endif
  }
  
  if (insideCode==kOutside) 
    belongVolume = false;

  else if (insideCode==kSurface) {
    belongVolume = false;
    reg = -100;
  }
  else {
    // Search downwards in daughter volumes
    // until deepest containing volume found
    //
    // 3 Cases:
    //
    // o Parameterised daughters
    //   =>Must be one G4PVParameterised daughter & voxels
    // o Positioned daughters & voxels
    // o Positioned daughters & no voxels
    
    belongVolume = true;    
    noResult=true;  
    
    do {
      // Determine `type' of current mother volume
      targetPhysical=fHistory->GetTopVolume();
      targetLogical=targetPhysical->GetLogicalVolume();
      switch(CharacteriseDaughters(targetLogical)) {
      case kNormal:
	if (targetLogical->GetVoxelHeader()) {
	  noResult=fvoxelNav.LevelLocate(*fHistory,
					 fBlockedPhysicalVolume,
					 fBlockedReplicaNo,
					 globalPoint,
					 pGlobalDirection,
					 fLocatedOnEdge,
					 localPoint);
	}
	else {
	  noResult=fnormalNav.LevelLocate(*fHistory,
					  fBlockedPhysicalVolume,
					  fBlockedReplicaNo,
					  globalPoint,
					  pGlobalDirection,
					  fLocatedOnEdge,
					  localPoint);
	}
	break;
      case kReplica:
	noResult=freplicaNav.LevelLocate(*fHistory,
					 fBlockedPhysicalVolume,
					 fBlockedReplicaNo,
					 globalPoint,
					 pGlobalDirection,
					 fLocatedOnEdge,
					 localPoint);
	break;
      case kParameterised:
	noResult=fparamNav.LevelLocate(*fHistory,
				       fBlockedPhysicalVolume,
				       fBlockedReplicaNo,
				       globalPoint,
				       pGlobalDirection,
				       fLocatedOnEdge,
				       localPoint);
	break;
      } //switch
      
      // LevelLocate search in the first daughter level. 
      // LevelLocate returns noResult=true if it finds a daughter volume
      // in which globalPoint is inside (or on the surface). So point
      // doesn't belong only to mother volume ==> belongVolume=false.
      
      if (noResult) {
	// The blocked volume no longer valid - it was for another level
	fBlockedPhysicalVolume= 0;
	fBlockedReplicaNo= -1;
	belongVolume=false;
      }
    } while (noResult);
    
    //on exit targetPhysical is the volume globalPoint belongs to;
    int volIndex = ~0;
    for (unsigned int i=0; i<pVolStore->size(); i++)
      if ((*pVolStore)[i] == targetPhysical) 
	volIndex = i;     
    if (volIndex==(~0)) {
      std::cerr << "FLUGG: Problem in routine WrapG1 tryingto find volume after step" << std::endl;
      exit(-999);
    }
    //reg = int(pVolStore->index(targetPhysical))+1;
    reg = volIndex+1;
    
  }
  
  delete fHistory;
  return belongVolume;
}

EVolume CharacteriseDaughters(const G4LogicalVolume *pLog)
{
  EVolume type;
  EAxis axis;
  int nReplicas;
  double width,offset;
  bool consuming;
  G4VPhysicalVolume *pVol;
  
  if (pLog->GetNoDaughters()==1) {
    pVol=pLog->GetDaughter(0);
    if (pVol->IsReplicated()) {
      pVol->GetReplicationData(axis,nReplicas,width,offset,consuming);
      type=(consuming) ? kReplica : kParameterised;
    }
    else
      type=kNormal;
  }
  else
    type=kNormal;
  return type;
}
*/

