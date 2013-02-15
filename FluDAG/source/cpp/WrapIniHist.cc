
// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapIniHist.hh - Sara Vanini
//
// Wrapper for reinitialization of FluggNavigator history.
//
// modified 14/I/99
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
///////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
// #include "FGeometryInit.hh"
// #include "NavHistWithCount.hh"
// #include "G4NavigationHistory.hh"
// #include "FluggNavigator.hh"
// #include "globals.hh"

void inihwr(int& intHist)
{
//flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "============= INIHWR ==============" << std::endl;    
  std::cout << "Ptr History=" <<intHist<< std::endl;
#endif 
/*
  if(intHist==-1) {
    std::cout << "ERROR! This history has been deleted!" << std::endl;
    return;
  }
  else {
    //get NavHistWithCount,G4NavigationHistory,FGeometryInit,
    //FluggNavigator pointers
    NavHistWithCount* ptrNavHistCount=
      reinterpret_cast<NavHistWithCount*>(intHist);
    G4NavigationHistory* ptrNavHist=ptrNavHistCount->GetNavHistPtr();
    static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
    FluggNavigator* ptrNavig = ptrGeoInit->getNavigatorForTracking();
    
    //reinitialize navigator history 
    ptrNavig->UpdateNavigatorHistory(ptrNavHist);
    
    //update utility histories: touch,temp, and reset old history
    ptrGeoInit->UpdateHistories(ptrNavHist,0);
    
    //save new history in jrLtGeant if not present
    int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
    int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
    
    bool intHistInJrLtGeant = false;
            bool oldversion = false;
    //	 bool oldversion = true;
    if ( oldversion ) {
    for(int h=0; h<=LttcFlagGeant; h++)
      if(jrLtGeant[h]==intHist) 
	intHistInJrLtGeant = true; }
    else if(LttcFlagGeant>=0) {
 
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "* CONHWR call to check and delete histories in jrLtGeant[]:"
	   << std::endl;
#endif 
    
    for(int ind=0; ind<=LttcFlagGeant; ind++) {
      int incrCount1=-1;
      conhwr(jrLtGeant[ind],&incrCount1);
      jrLtGeant[ind]=-1;
          } 
    LttcFlagGeant=-1;
      } 
    if(!intHistInJrLtGeant) {  
      LttcFlagGeant += 1; 
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "* WrapIniHist:  LttcFlagGeant increase" << LttcFlagGeant <<std::endl;
#endif
      ptrGeoInit->SetLttcFlagGeant(LttcFlagGeant);
      
      jrLtGeant[LttcFlagGeant]=intHist;
      
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "* CONHWR call to increment counter" << std::endl;
#endif 
      int incrCount=1;
      conhwr(jrLtGeant[LttcFlagGeant],&incrCount);
    }    
    
    //print history....
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "History reinitialized in:" << std::endl;
    std::cout << *ptrNavHist << std::endl;
    ptrGeoInit->PrintJrLtGeant();
#endif
  }
*/
}





