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
#define DEBUG 1

extern "C" {
 void flukam(const int &GeoFlag);
}

// Perform the task of argument and error checking
void checkArgs(std::string infile, int numargs, char* argv[], bool& flukarun);
// char *prefixFilename(char *cfile, bool running_with_fluka);
std::string prefixFilename(char *cfile, bool running_with_fluka);

int main(int argc, char* argv[]) 
{
  bool flukarun = false;
  
  // Default h5m filename is for fluka runs
  std::string infile = "test.h5m";
  std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
  checkArgs(infile, argc, argv, flukarun);

  char *fileptr;
  fileptr = &infile[0];
  // std::strcpy(fileptr, infile.c_str());
  // char *locatedFile = prefixFilename(fileptr, flukarun);
  std::string locatedFile = prefixFilename(fileptr, flukarun);
  std::cerr << "prefixed file is " << locatedFile << std::endl;

  bool success  = checkFile(locatedFile);
  if (!success)
  {
     exit(EXIT_FAILURE);
  }


  // Load the h5m file, init the obb tree;  flukarun changes the expected
  // relative location of the file
  int max_pbl = 1;
  char *ptr = new char[locatedFile.length()+1];
  std::strcpy(ptr, locatedFile.c_str());
  cpp_dagmcinit(ptr, 0, max_pbl); 
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
  std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
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
      std::cerr << "Writing materials to " << infile << std::endl;
      flukarun = false;
  }
} 


//---------------------------------------------------------------------------//
// prefixFilename(..)
//---------------------------------------------------------------------------//
// If running_with_fluka, prefix "../" to the filename, as fluka runs from a subdirectory
//char *prefixFilename(char *cfile, bool running_with_fluka)
std::string prefixFilename(char *cfile, bool running_with_fluka)
{
  std::string prefixedFilename; 
  // Prefix
  if (running_with_fluka)  // h5m file is one level up
  {
     prefixedFilename="../";
  }
  else // file is in same directory as executable
  {
     prefixedFilename="";
  }
  prefixedFilename.append(cfile);
  if (DEBUG)
  {
  	std::cerr << "\nmy file is " << prefixedFilename << "\n" << std::endl;
  }
  char *myfile;
  myfile = &prefixedFilename[0];
  // return myfile;
  return prefixedFilename;
}  



