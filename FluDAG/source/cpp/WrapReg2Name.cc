
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
  std::cout << "============= RG2NWR ==============" << std::endl;    
  std::cout << "mreg=" << mreg << std::endl;
  char * vvname;
  region2name(mreg, vvname);
  Vname = vvname;
/*
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();

  G4String VVname=  ptrGeoInit->fRegionNameMap[mreg];

  VVname.resize(8,' ');
  for (int i=0; i<8; i++){
    Vname[i] = VVname.at(i);
    //    std::cout <<i<< " Vname(i) " << Vname[i] << std::endl; }
  }
*/
  std::cout << "reg2nmwr: Vname " << Vname<< std::endl;  
}





