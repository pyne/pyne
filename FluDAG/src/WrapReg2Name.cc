
// FluDAG tag 

///////////////////////////////////////////////////////////////////
//
// WrapReg2Name.cc -  Paola Sala
//
// Wrapper for getting region name corresponding to given region number
//
// created 8-jul-2009
//
///////////////////////////////////////////////////////////////////

#include "DagWrappers.hh"
#include "DagWrapUtils.hh"


void rg2nwr(const int& mreg, const char* Vname)
{
  std::cerr << "============= RG2NWR ==============" << std::endl;    
  std::cerr << "mreg=" << mreg << std::endl;
  char * vvname;
  region2name(mreg, vvname);
  Vname = vvname;
  std::cerr << "reg2nmwr: Vname " << Vname<< std::endl;  
  return;
}





