
// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapSavHist.hh - Sara Vanini
//
// Wrapper for saving current navigation history (fCheck=default) 
// and returning its pointer. If fCheck=-1 copy of history pointed 
// by intHist is made in NavHistWithCount object, and its pointer 
// is returned. fCheck=1 and fCheck=2 cases are only in debugging 
// version: an array is created by means of FGeometryInit functions
// (but could be a static int * ptrArray = new int[10000] with 
// file scope as well) that stores a flag for deleted/undeleted 
// histories and at the end of event is checked to verify that 
// all saved history objects have been deleted.
//
// modified 6/III/99: history check array implemented
// modified 14/IV/00: fCheck=-1 case modified
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
//////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
// #include "FGeometryInit.hh"
// #include "NavHistWithCount.hh"
// #include "G4TouchableHistory.hh"
// #include "globals.hh"


int isvhwr(const int& fCheck, const int& intHist)
{
  //flag
#ifdef DAGGEOMETRY_DEBUG
  std::cout << "============= ISVHWR ==============" << std::endl;    
  std::cout << "fCheck=" << fCheck << std::endl;
  if(fCheck==-1) 
    std::cout << "intHist=" << intHist  << std::endl;
#endif 
 /* 
  //Geoinit pointer
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  static int j=0;
  
  switch (fCheck) {
  case 1:
    {
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "Start event." << std::endl;
#endif
      
      return 0;
    }
    
  case 2:
    {
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "End event. Check histories in debug-array..." << std::endl;
#endif
*/      
      //check that all fluka-banked histories have been delated
      // commented by A.Solodkov

      /*
	int * ptrArray = ptrGeoInit->GetHistArray();
	
	for(int k=0;k<j;k++) {
	NavHistWithCount* ptrNavHistCount=reinterpret_cast
	<NavHistWithCount*>(ptrArray[k]);
	
	if(ptrArray[k] && !ptrNavHistCount->GetDelateFlag()) 
	std::cout << "WARNING! History pointed by " <<ptrArray[k]<<
	" has not been deleted at end of event" << std::endl;
	}
	
	//reinitialise debug histories array
	for(int i=0;i<1000000;i++) ptrArray[i]=0;
      */
/*
      j=0;
      
      return 0;
    }
    
  case -1:
    {
      //get history object from intHist
      NavHistWithCount* ptrNavHistCountCopy = 
	reinterpret_cast<NavHistWithCount*>(intHist);
      G4NavigationHistory* ptrNavHistCopy = 
	ptrNavHistCountCopy->GetNavHistPtr();
      
      //copy history in another NavHistWithCount object
      //and save index of check-array
      NavHistWithCount * ptrNavHistCount = 
	new NavHistWithCount(*(ptrNavHistCopy));
      ptrNavHistCount->SaveCheckInd(j);
      int intHistCopy = long(ptrNavHistCount);
      
      //store ptr in array
      // commented by PS
      //int * ptrArray = ptrGeoInit->GetHistArray();
      //ptrArray[j]=intHistCopy;
      j+=1;
      
#ifdef DAGGEOMETRY_DEBUG
      std::cout << "Copying history..." << std::endl;
      std::cout << "Ptr History-copy =" <<intHistCopy<< std::endl;
      std::cout<<*ptrNavHistCopy<< std::endl;
#endif 
      
      return intHistCopy;
    }
    
  default:
    {
      //save G4NavigationHistory and index of check-array
      NavHistWithCount *ptrNavHistCount = 
	new NavHistWithCount(*(ptrGeoInit->GetTempNavHist()->GetHistory()));
      
      ptrNavHistCount->SaveCheckInd(j);
      int histInt = long(ptrNavHistCount);
      
      //store ptr in array
      // comented by PS
      //int * ptrArray = ptrGeoInit->GetHistArray();
      // ptrArray[j]=histInt;
      j+=1;
      
#ifdef DAGGEOMETRY_DEBUG
      //TouchableHistory 
      G4TouchableHistory * ptrTouchHist = ptrGeoInit->GetTouchableHistory();
      std::cout << "Saving history..." << std::endl;
      std::cout << "Ptr saved History=" << histInt << std::endl;
      std::cout << *(ptrTouchHist->GetHistory()) << std::endl;
#endif 
      
      return histInt;
    }
  }
*/
}






