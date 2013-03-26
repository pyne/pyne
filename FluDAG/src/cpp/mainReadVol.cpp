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
#include "fluka_funcs.h"
#include "dagmc_utils.hpp"

ErrorCode readVol(char *fileptr, const std::string* override_output_filename=NULL );

int main(int argc, char* argv[]) 
{
  ErrorCode code;
  if (argc < 2)
  {
     std::cerr << "Usage:  " << argv[0] << " filename.h5m" << std::endl;
     exit(0);
  }
 
  int max_pbl = 1;
  bool flukarun = false;
  char *fileptr = argv[1];
  cpp_dagmcinit(fileptr, 0, max_pbl,flukarun); 

  code = readVol(argv[1]);
  return 0;
}
