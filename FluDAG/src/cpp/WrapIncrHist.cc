
// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapIncrHist.hh - Sara Vanini
//
// Wrapper for updating  secondary particles history counter. 
// If counter=0 the history is deleted. 
//
//////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

void conhwr(int& intHist, int* incrCount)
{
  std::cerr << "============= CONHWR ==============" << std::endl;    
  std::cerr << "Ptr History = " << intHist << std::endl;
  incrCount++;
  std::cerr << "Counter = " << incrCount << std::endl;
  std::cerr << "============= Out of CONHWR ==============" << std::endl;    
  
  return;
 }










