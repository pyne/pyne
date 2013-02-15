#include "fluka_funcs.h"

#include "moab/Interface.hpp"
#include "DagMC.hpp"
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string.h>


using namespace moab;

#define flukam flukam_
#define DAG DagMC::instance()

extern "C" {
 void flukam(const int &GeoFlag);
}


int main(int argc, char* argv[]) {

  bool flukarun = false;

  // std::string infile = "model_complete.h5m";
  std::string infile = "test.h5m";

  char * fileptr = new char [infile.length()+1];
  std::strcpy(fileptr, infile.c_str());
  // No filename => do a fluka run using test.h5m in higher directory
  if (argc < 2) {
                // Tell the user how to run the program
                std::cerr << "Using " << infile << std::endl;
                std::cerr << "   or call: " << argv[0] << " h5mfile" << std::endl;
                /* "Usage messages" are a conventional way of telling the user
                 * how to run a program if they enter the command incorrectly.
                 */
                flukarun = true;
  }
  else  // Give a file name to write out the material file and stop
  {
      std::strcpy(fileptr, argv[1]);
      std::cerr << "Using " << fileptr << std::endl;
      flukarun = false;
  }
  int max_pbl = 1;
  // Load the h5m file, init the obb tree   
  cpp_dagmcinit(fileptr, 0, max_pbl, flukarun);
  
  if (!flukarun)
  {
    std::string lcad = "mat.inp";
    fludagwrite_assignma(lcad);
  }

//flag for geometry:
// 1 for GEANT4
// 0 for FLUKA
// 2 for Rubia
// 3 for Dagmc ?
    const int flag = 1;

//call fortran
// Temporarily comment out while testing the writing of FlukaMat
    if (flukarun)
    {
       flukam(flag);
    }

//end
  return 0;
}




