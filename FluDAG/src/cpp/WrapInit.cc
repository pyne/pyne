
///////////////////////////////////////////////////////////////////
//
// WrapInit.hh - 
//
// Wrapper for geometry initialisation.
//
//////////////////////////////////////////////////////////////////


#include "DagWrapUtils.hh"
#include "DagWrappers.hh"

using namespace moab;

void jomiwr(int & nge, const int& lin, const int& lou, int& flukaReg)
{
  std::cerr << "================== JOMIWR =================" << std::endl;
  // return FLUGG code to fluka
  // nge = 3;

  //Original comment:  returns number of volumes + 1
  unsigned int numVol = DAG->num_entities(3);
  flukaReg = numVol;
	
  std::cerr << "Number of volumes: " << flukaReg << std::endl;
  std::cerr << "================== Out of JOMIWR =================" << std::endl;

  return;
}
