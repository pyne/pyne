

///////////////////////////////////////////////////////////////////
//
// WrapG1.hh - Sara Vanini
//
// Wrapper for geometry tracking: returns approved step of 
// particle and alla variables that fluka G1 computes.
//
// modified 11/III/99 : included jrLtGeant array for storing 
// lattice histories; fixed jump-on-boundaries and back-scattering
// modified 20/IV/99 : history counters of jrLt array are
// incremented when created and decremented when check is required:
// this for not deleting histories twice.
// modified 18/X/99 : LocateGlobalPointWithinVolume used when position
// is changed from last located point.
// modified 1/III/00 : update utilities histories when crossing 
// identical volume boundaries.
//
// modified 22/III/00 : fixed LttcFlag and jrLt return values.
// modified 12/VI/00 : end-step on Boundary bug fixed.
// modified 5/VII/00 : boundary not seen by G4 geometry bug fixed.
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
/////////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"
// #include "FGeometryInit.hh"
// #include "NavHistWithCount.hh"
// #include "G4VPhysicalVolume.hh"
// #include "G4NavigationHistory.hh"
// #include "G4GeometryTolerance.hh"
// #include "FluggNavigator.hh"
// #include "G4PhysicalVolumeStore.hh"
// #include "globals.hh"


void g1wr(double& pSx, double& pSy, double& pSz, double* pV,
          int& oldReg, const int& oldLttc, double& propStep,
          int& nascFlag, double& retStep, int& newReg,
          double& saf, int& newLttc, int& LttcFlag,
          double* sLt, int* jrLt)
{
  //flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout<<"============= G1WR =============="<<std::endl;    
#endif 
  
  //static int count=0;
  //count+=1;
  //std::cout<<"contatore G1="<<count<<std::endl;
  
  ///////////////////////// INPUT ///////////////////////////
  //Geoinit, Navigator, TouchableHistory pointers
/*
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  FluggNavigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
  G4TouchableHistory * ptrTouchHist = ptrGeoInit->GetTouchableHistory();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
*/
  
  //setting variables (and dimension: Fluka uses cm.!)
/*
  G4ThreeVector partLoc(pSx,pSy,pSz);
  partLoc *= 10.0; // in millimeters!
  static G4ThreeVector partLocOld = partLoc;
  static G4ThreeVector oldLocalPoint = 
    ptrNavig->GetGlobalToLocalTransform().TransformPoint(partLocOld);
  
  G4ThreeVector pVec(pV[0],pV[1],pV[2]);
  const double physStep=double(propStep*10.);
  
  int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
 */ 
  
#ifdef DAGGEOMETRY_DEBUG
/*
  std::cout.precision(10);
  std::cout << "Position (cm):" << pSx << "," << pSy << "," << pSz << std::endl;
  std::cout << "Direction: "    << pVec << std::endl;
  std::cout << "Proposed step :"<< propStep << std::endl;
*/
#endif 



  ///////////////////// FLUKA/G4 REGIONS COMPARISON //////////////////////
  //get oldReg pointer and G4volume of previous step end-point pointer
/*
  G4VPhysicalVolume * ptrOldReg = (*pVolStore)[oldReg-1];
  G4VPhysicalVolume * lastVolume = ptrTouchHist->GetVolume();
  int touVolCpNb = 0;
  if (lastVolume) 
    touVolCpNb = lastVolume->GetCopyNo();
 */ 
#ifdef  DAGGEOMETRY_DEBUG
/*
  int fluVolCpNb = ptrOldReg->GetCopyNo();
  std::cout << "Fluka volume before step: " << ptrOldReg->GetName()
	 << "," << fluVolCpNb<< std::endl;
  if(lastVolume) 
    std::cout << "G4 Touch Hist volume: "
	   << lastVolume->GetName() << "," << touVolCpNb << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
*/
#endif 


  //if volume is changed, this is a new particle tracking, or fluka tries
  //to reach a boundary more softly with the same particle. In this case
  //fluka restart tracking from old history, in general. For tracking in 
  //lattice-volume fluka goes back to one of the previous lattice volumes. 
  //Check if ptrOldReg is equal to old history, or to one of the N lattice 
  //volumes stored in jrLt. Then reinitialise step-histories! Otherwise 
  //must relocate.
  //NB. jrLtGeant stores lattice volume histories until LttcFlag==-1, 
  //then all histories are checked and deleted, and array is reinitialised 
  
  
/*
  int haveHistNb = -1;
  int newRegErr=0;
  int indHist = LttcFlagGeant;
  int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
  
  
  G4NavigationHistory* ptrLttcHist=0;
  if(oldLttc) 
    ptrLttcHist = reinterpret_cast<NavHistWithCount*>(oldLttc)->GetNavHistPtr();
  
  
  while(indHist>=0 && haveHistNb==-1) {
#ifdef  DAGGEOMETRY_DEBUG
    std::cout<<"Searching in jrLtc...."<<std::endl;
#endif 
    if(oldLttc==jrLtGeant[indHist]) 
      haveHistNb=indHist;
    indHist-=1;
  }
  
*/
/*
  if(haveHistNb!=-1) {
    //fluka history found in jrLtGeant
    if(haveHistNb<LttcFlagGeant) {
#ifdef  DAGGEOMETRY_DEBUG
      std::cout<<"* Fluka reaches boundary more softly..."<<std::endl;
      std::cout<<"* Re-initializing FluggNavigator history"<<std::endl;
      std::cout<<"* and updating step-histories"<<std::endl;
#endif 
      
      ptrNavig->UpdateNavigatorHistory(ptrLttcHist);
      ptrGeoInit->UpdateHistories(ptrLttcHist,0);
    }
#ifdef  DAGGEOMETRY_DEBUG
    if(haveHistNb==LttcFlagGeant) 
      std::cout << "Continuing step...." << std::endl;
#endif 
    jrLt[0]=oldLttc;
  }
  else {
    //not found fluka history in jrLttGeant!
    std::cout << "* ERROR! Geometry not correctly initialised in fluka history!"
	   << std::endl; 
    std::cout << "Position (cm):" << pSx << "," << pSy << "," << pSz << std::endl;
    std::cout << "Direction: "    << pVec << std::endl;
    std::cout << "Proposed step :"<< propStep << std::endl;
    int fluVolCpNb = ptrOldReg->GetCopyNo();
    std::cout << "Fluka volume before step: " << oldReg << " " << ptrOldReg->GetName()
	   << "," << fluVolCpNb<< std::endl;
    if(lastVolume) 
      std::cout << "G4 Touch Hist volume: "
	     << lastVolume->GetName() << "," << touVolCpNb << std::endl;

    //relocation!
    ptrNavig->LocateGlobalPointAndUpdateTouchable(partLoc,pVec,ptrTouchHist,true);
    
    std::cout << "* ATTENTION: point relocation in: "
	   << ptrTouchHist->GetVolume()->GetName() << std::endl;
    
    ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),1);
    
    
    //save new history in jrLt[0] and increment its counter 
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* ISVHWR call to store new NavHistWithCount in jrLt[0]"
	   << std::endl;
#endif 
    jrLt[0]=isvhwr(0,0);

#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* CONHWR call to increment counter" << std::endl;
#endif 
    int incrCount2=1;
    conhwr(jrLt[0],&incrCount2);
  }
  
    if(ptrTouchHist->GetVolume() != ptrOldReg) {
      std::cout << "* Warning! Point not in fluka volume!" << std::endl;
    std::cout << "Position (cm):" << pSx << "," << pSy << "," << pSz << std::endl;
    std::cout << "Direction: "    << pVec << std::endl;
    std::cout << "Proposed step :"<< propStep << std::endl;
    newRegErr=-3;
    }
  
  //jrLtGeant - history check: decrement counter and delete 
  //histories, if LttcFlag=-1, then reinitialise array with -1. 
  if(LttcFlag==-1 && LttcFlagGeant>=0) {
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* CONHWR call to check and delete histories in jrLtGeant[]:"
	   << std::endl;
#endif 
    
    for(int ind=0; ind<=LttcFlagGeant; ind++) {
      int incrCount1=-1;
      if(jrLtGeant[ind]!=jrLt[0]) conhwr(jrLtGeant[ind],&incrCount1);
      //      if(cpLtGeant[ind]!=-1) conhwr(cpLtGeant[ind],&incrCount1);
      jrLtGeant[ind]=-1;
      //      cpLtGeant[ind]=-1;
    } 
    LttcFlagGeant=-1;
  }
  
  
  //update jrLt and sLt arrays   
  int end = 99;
  if(LttcFlag>=0) 
    end=LttcFlag;

  // Added by A.Solodkov
  if (end>=100) {
    std::cout << "Problems in WrapG1 routine" << std::endl;
    std::cout << "Index LttcFlag=" << end << " is outside array bounds" << std::endl;
    std::cout << "Better to stop immediately !" << std::endl;
    exit(1);
  }
  //jrLt re-initialization with -1 (jrLt[0] is already set)
  for(int vv=1;vv<=end;vv++) 
    jrLt[vv]=-1;
  //sLt re-initialization
  for(int vs=0;vs<=end;vs++) 
    sLt[vs]=0;
  
  LttcFlag=0;
  
  
  /////////////////////  COMPUTE STEP  ////////////////////////
  //update Navigator private flags and voxel stack if point is 
  //changed from last located point (otherwise troubles come 
  //when fluka changes location or particle because G4 computes 
  //from last located point).
  G4ThreeVector newLocalPoint = ptrNavig->GetGlobalToLocalTransform().TransformPoint(partLoc);
  double moveLenSq = (newLocalPoint-oldLocalPoint).mag2();
  // new in G4.9.0: kCarTolerance does not exist anymore
  // one needs to get the tolerance from G4GeometryTolerance class
  double carTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    
  if(moveLenSq>=carTolerance*carTolerance) 
    ptrNavig->LocateGlobalPointWithinVolume(partLoc);
  
  
  //compute step and new location
  newReg = oldReg;
  double appStep = 0;
  double safety = 0;
  G4bool onBoundaryRound = false;
  G4bool crossBound = false;
  double physStepTmp = physStep;
  G4bool flagError = false;
  
  while( (oldReg==newReg && appStep<physStepTmp) || onBoundaryRound ) {
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* Computing step..." << std::endl;
    std::cout <<"  appStep "<< appStep << " physStepTmp  "<< physStepTmp  << std::endl;
    std::cout <<" onBoundaryRound  " <<  onBoundaryRound << std::endl;
#endif
    
    //update variables
    oldReg = newReg;
    
    if(onBoundaryRound) {
      physStepTmp=10.e-10;
      //compute step and location: returns newReg
      newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
			     safety,onBoundaryRound,flagError,oldReg);
      if (onBoundaryRound && oldReg== newReg) {
	// again! means the navigator gives strange answer, stop here to avoid infinite loop
	std::cout << " onBoundaryRound loop  " << oldReg <<std::endl;
	onBoundaryRound=false;
	  physStepTmp=0.;
      }
      physStepTmp=0.;
      crossBound=true;
    }
    else {
      physStepTmp -= appStep;
      //compute step and location: returns newReg
      newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
			     safety,onBoundaryRound,flagError,oldReg);
    }
    
    G4bool EqualHist;
    EqualHist = EqualHistories(ptrTouchHist->GetHistory(),
			       ptrGeoInit->GetTempNavHist()->GetHistory());
    
    if(!EqualHist && flagError) {
      pVec=-pVec;
      newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
                             safety,onBoundaryRound,flagError,oldReg); 
      pVec=-pVec;
      physStepTmp+=1;
      newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
                             safety,onBoundaryRound,flagError,oldReg);
      
      EqualHist = EqualHistories(ptrTouchHist->GetHistory(),
				 ptrGeoInit->GetTempNavHist()->GetHistory()); 
    }
    
    //update sLt
    double pas = double(appStep);
    sLt[LttcFlag] += pas/cm;
    safety = (oldReg!=newReg)?0.:safety; 
    
    if(!EqualHist) {
      //end-step point is on boundary between volumes;
      //==> update step-histories, save new NavHistWithCount 
      //and save its pointer in jrLt; save step in sLt.
      
      //set onBoundaryRound=false to avoid re-compute step!
      onBoundaryRound=false;
      
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "* History is changed!" << std::endl;
      std::cout << "* updating step-histories, jrLt, LttcFlag" << std::endl;
#endif
      
      ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),2);
      
      LttcFlag += 1;
      
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "* ISVHWR call to store new NavHistWithCount in jrLt" 
	     << std::endl;
#endif 
      
      jrLt[LttcFlag] = isvhwr(0,0);
      
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "* CONHWR call to increment counter" << std::endl;
#endif 
      int incrCount3=1;
      conhwr(jrLt[LttcFlag],&incrCount3);
      
      sLt[LttcFlag] = sLt[LttcFlag-1];
    }
  }
  
  //////////////////////   OUTPUT   //////////////////////////
  //If back-scattering occured, and fluka is in the wrong region, return -3.
  //(N. B. Boundary between replicans are not seen when physStep=distance 
  //particle-boundary: in this case step=kInfinity and history is unchanged. 
  //Following step is =0 and then history changes.) 
  if(nascFlag<0 && !appStep && physStep && newReg!=oldReg && !crossBound) {
    //don't need to compare histories because boundary between
    //identical volumes in different replicans are not seen
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* Back-scattering!" << std::endl;
    std::cout << " nascFlag  " << nascFlag << " appStep " << appStep  << " physStep  " <<   physStep << std::endl;
    std::cout  << " newReg  " << newReg  << "oldReg  " << oldReg  << "crossBound  " << crossBound  << std::endl;
#endif
    // PS commented out on 12-4-2011     
    //    newReg=-3;
  }
  else if(newRegErr<0) 
    {
      std::cout << " new Reg err " <<   newReg << std::endl;
      newReg=newRegErr;
    }
  
  //compute output variables (in cm.!)
  //final step
  retStep = sLt[LttcFlag];
  
  //safety (Fluka sottracts a bit to safety to be sure 
  //not to jump on a boundary)
  double s = double(safety);
  s -= s*3.0e-09;
  saf=s/cm; 
  
  //update wrapper utility variables
  //copy jrLt in jrLtGeant
  int start=0;
  if(haveHistNb!=-1 && LttcFlagGeant!=-1) 
    start=1;
  for(int lt=start;lt<=LttcFlag;lt++)
    jrLtGeant[LttcFlagGeant+1+lt-start]=jrLt[lt];
  LttcFlagGeant+=(1+LttcFlag-start);
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "* WrapG1:  LttcFlagGeant" << LttcFlagGeant <<std::endl;
#endif
  newLttc = jrLt[LttcFlag];
  ptrGeoInit->SetLttcFlagGeant(LttcFlagGeant);
  
  partLocOld=partLoc;
  oldLocalPoint = ptrNavig->GetGlobalToLocalTransform().TransformPoint(partLocOld);
  
  //compute new position
  G4ThreeVector oldPos = G4ThreeVector(pSx,pSy,pSz);
  G4ThreeVector newPos = oldPos + retStep*pVec;
  
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "New position: " << newPos << std::endl;
  std::cout << "Output region: " << newReg << std::endl;
  std::cout << "G4 safety (cm): " << (safety*0.1) << std::endl;
  std::cout << "Fluka safety (cm): " << saf << std::endl;
  std::cout << "Approved step: " << retStep << std::endl;
  std::cout << "LttcFlag = " << LttcFlag << std::endl;
  for(int i=0;i<=LttcFlag+1;i++) {
    std::cout << "jrLt[" << i <<"]=" << jrLt[i] << std::endl;
    std::cout << "sLt[" << i <<"]=" << sLt[i] << std::endl;
  }
  
  std::cout << "LttcFlagGeant =" << LttcFlagGeant << std::endl;
  for(int ib=0;ib<=LttcFlagGeant+1;ib++) {
    std::cout << "jrLtGeant[" << ib <<"]=" << jrLtGeant[ib] << std::endl;
  }
  std::cout << "newLttc= " << newLttc << std::endl;
#endif
*/
}




