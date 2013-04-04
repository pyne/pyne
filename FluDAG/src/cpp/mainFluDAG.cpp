#include "fludag_utils.h"
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

// Perform the task of argument and error checking
void checkArgs(std::string infile, int numargs, char* argv[], bool& flukarun);

int main(int argc, char* argv[]) 
{
  bool flukarun = false;
  
  // Default h5m filename is for fluka runs
  std::string infile = "test.h5m";
  checkArgs(infile, argc, argv, flukarun);

  char *myfile = checkInput(infile, flukarun);

  // Load the h5m file, init the obb tree;  flukarun changes the expected
  // relative location of the file
  int max_pbl = 1;
  cpp_dagmcinit(myfile, 0, max_pbl); 
  if (!flukarun)
  {
    std::string lcad = "mat.inp";
    fludagwrite_assignma(lcad);
    fludagwrite_mat("mat1.inp");
  }
  else // call flukarun
  {
//flag for geometry:
// 1 for GEANT4
// 0 for FLUKA
// 2 for Rubia
// 3 for Dagmc ?

       const int flag = 1;
       flukam(flag);
    }

  return 0;
}

// Perform argument and argument-dependent checks
void checkArgs(std::string infile, int numargs, char* argv[], bool& flukarun)
{
  // No filename => do a fluka run using test.h5m in higher directory
  if (numargs < 2) 
  {
     // Tell the user how to run the program
     std::cerr << "Using " << infile << std::endl;
     std::cerr << "   or call: " << argv[0] << " h5mfile" << std::endl;
     // "Usage messages" are a conventional way of telling the user
     // how to run a program if they enter the command incorrectly.
     flukarun = true;
  }
  else  // Given a file name, write out the material file and stop
  {
      std::cerr << "Using " << infile << std::endl;
      flukarun = false;
  }
} 



