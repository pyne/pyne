#include "fluka_funcs.h"

#include "moab/Interface.hpp"
#include "DagMC.hpp"
#include <iostream>
#include <stdlib.h>
#include <string.h>


using namespace moab;

#define flukam flukam_
#define DAG DagMC::instance()

extern "C" {
 void flukam(const int &GeoFlag);
}

 // ToDo: move to header
 ErrorCode readVol(char *file);

int main() {

  int max_pbl = 1;
  // ToDo: Hardcode the filename temporarily
  std::string infile = "test.h5m";
  char *fileptr;
  fileptr = &(infile[0]);
  
  // Load the h5m file, init the obb tree   
  cpp_dagmcinit(fileptr, 0, max_pbl);
  
  ErrorCode returncode;
  // ToDo: fileptr does not need to be passed, just used for a print statement
  returncode = readVol(fileptr);
  std::vector<std::string> detected;
  DAG->detect_available_props(detected);
  DAG->parse_properties (detected);
  std::cout << detected.size() << " metadata properties detected:" << std::endl;
  for( std::vector<std::string>::iterator kwdp = detected.begin();
          kwdp != detected.end(); ++kwdp )
  {
       std::string keyword = *kwdp;
       std::cout << "    " << keyword << std::endl;
  }
// Answer for test.h5m:
///////////////////////

// 5 metadata properties detected:
//    graveyard
//    mat
//    rho
//    surf.flux
//    tally
/////////////////

//flag for geometry:
// 1 for GEANT4
// 0 for FLUKA
// 2 for Rubia
// 3 for Dagmc ?
    const int flag = 1;

//call fortran
// Temporarily comment out while testing the writing of FlukaMat
    // flukam(flag);

//end
  return 0;
}




