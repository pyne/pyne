//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/mainReadVol.cpp
 * \author Julie Zachman 
 * \date   Fri Mar 8 2013 
 * \brief  Read properties from an h5m file
 * \note   
 */
//---------------------------------------------------------------------------//
// $Id: 
//---------------------------------------------------------------------------//
#include "fludag_utils.h"
#include "fluka_funcs.h"
#include "dagmc_utils.hpp"

ErrorCode readVol(char *fileptr, const std::string* override_output_filename=NULL );

// Perform the task of argument and error checking
void checkArgs(int numargs, char* argv[]);

int main(int argc, char* argv[]) 
{
  ErrorCode code;
  checkArgs(argc, argv);
 
  int max_pbl = 1;
  char *fileptr = argv[1];
  cpp_dagmcinit(fileptr, 0, max_pbl); 

  code = readVol(argv[1]);
  return 0;
}

// Arg check
void checkArgs(int numargs, char* args[])
{
  if (numargs < 2)
  {
     std::cerr << "Usage:  " << args[0] << " filename.h5m" << std::endl;
     exit(EXIT_FAILURE);
  }
}
