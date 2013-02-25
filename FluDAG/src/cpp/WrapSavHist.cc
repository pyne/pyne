
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
//
//////////////////////////////////////////////////////////////////


#include "DagWrappers.hh"
#include "DagWrapUtils.hh"


int isvhwr(const int& fCheck, const int& intHist)
{
  std::cerr << "============= ISVHWR ==============" << std::endl;    
  std::cerr << "fCheck=" << fCheck << std::endl;
  if(fCheck==-1) 
    {
      std::cerr << "intHist=" << intHist  << std::endl;
    }

  return 1;
}

