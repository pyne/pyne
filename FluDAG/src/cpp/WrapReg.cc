
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
//
///////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"

void rgrpwr(const int& flukaReg, const int& ptrLttc, int& g4Reg,
            int* indMother, int* repMother, int& depthFluka)
{
  std::cerr << "============= RGRPWR ==============" << std::endl;    
  std::cerr << "ptrLttc=" << ptrLttc << std::endl;
  return;
}





