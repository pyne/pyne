
// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapIncrHist.hh - Sara Vanini
//
// Wrapper for updating  secondary particles history counter. 
// If counter=0 the history is deleted. 
//
// modified 14/I/9999
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
//////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
// #include "FGeometryInit.hh"
// #include "NavHistWithCount.hh"
// #include "globals.hh"


void conhwr(int& intHist, int* incrCount)
{
//flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "============= CONHWR ==============" << std::endl;    
  std::cout << "Ptr History = " << intHist << std::endl;
#endif 
  
/*
  //get NavHistWithCount pointer
  if(intHist!=-1) {
    NavHistWithCount* ptrNavHistCount=
      reinterpret_cast<NavHistWithCount*>(intHist); 
 
    //for debugging...
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "Secondary counter=" << ptrNavHistCount->GetCount();
    if(*incrCount>0) 
      std::cout << "+" << *incrCount << std::endl;
    if(*incrCount<0) 
      std::cout<< *incrCount << std::endl;   
    if(*incrCount==0) 
      std::cout << std::endl; 
#endif 
    
    //update secondary particles counter
    ptrNavHistCount->UpdateCount(*incrCount);
    
    //delete history if counter=0 or if counter=-1
    int counter = ptrNavHistCount->GetCount();
#ifdef DAGGEOMETRY_DEBUG
    std::cout << "Counter = " << counter << std::endl;
#endif 
    if(!counter || counter==-1) {
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "Deleting Nav Hist object..." << std::endl;
#endif 
*/
      /*
	//for history checking....
	int index = ptrNavHistCount->GetCheckInd();
	static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
	int * ptrArray = ptrGeoInit->GetHistArray();
	ptrArray[index]=0; 
      */	  
 /*     
      delete ptrNavHistCount;
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "end delete" << std::endl;
#endif 
      intHist=-1;
    }
  }
*/
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "============= Out of CONHWR ==============" << std::endl;    
#endif 
 }










