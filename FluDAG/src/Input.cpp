//----------------------------------*-C++, Fortran-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/Input.cpp
 * \author Julie Zachman 
 * \date   Apr 5 2013 
 * \brief  Functions called by fluka
 * \note   Unittested
 */
//---------------------------------------------------------------------------//
#include "fludag_utils.h"
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string.h>

//---------------------------------------------------------------------------//
// checkFile(..)
//---------------------------------------------------------------------------//
/// Call inline function to check basic i/o
bool checkFile(std::string filename)
{
  return isFileReadable(filename) && isFileNonEmpty(filename);
}


//---------------------------------------------------------------------------//
// isFileReadable(..)
//---------------------------------------------------------------------------//
// Test the non-prefixed file because the testing is done prior to calling fluka
bool isFileReadable(std::string filename)
{
  
  char *fileptr = new char[filename.length() + 1];
  std::strcpy(fileptr, filename.c_str());

  return file_can_be_read(fileptr);
/*
  bool success =   std::ifstream(fileptr); 
  if (success)
  {
    std::cerr << "Input file: " << fileptr << " can be read." << std::endl;
  }
  return success;
*/
}

//---------------------------------------------------------------------------//
// isFileNonEmpty(..)
//---------------------------------------------------------------------------//
// Test the prefixed file because the testing is done prior to calling fluka
bool isFileNonEmpty(std::string filename)
{
  char *fileptr = new char[filename.length() + 1];
  std::strcpy(fileptr, filename.c_str());

  std::ifstream file(fileptr, std::ios::binary);
        
  file.seekg(0, std::ios::end);
  if (file.tellg() == 0)
  {
        return false; // file is empty
  }
  else
  {
        return true;
  }
}
