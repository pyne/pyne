
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
  std::cout << "================== JOMIWR =================" << std::endl;
  // return FLUGG code to fluka
  nge = 3;

 // std::cout << "\t *==> JOMIWR: Setting DAGMC::instance()..." << std::endl;
  
  
  //create Fluka material cards in flukaMat.inp file
  // std::cout << "\t *==> JOMIWR: Init fluka materials...moved to main jcz" << std::endl;
  //createFlukaMatFile();
  
  //Original comment:  returns number of volumes + 1
  unsigned int numVol = DAG->num_entities(3);
  flukaReg = numVol;
	
  std::cout << "Number of volumes: " << flukaReg << std::endl;
  std::cout << "================== Out of JOMIWR =================" << std::endl;
}
