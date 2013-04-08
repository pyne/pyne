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

#define DEBUG 1

//---------------------------------------------------------------------------//
// checkInput(..)
//---------------------------------------------------------------------------//
/// Top level input checker
// Precondition: argument number has already been checked
// Handle error conditions gracefully
// bool checkInput(char *locatedFile)
// {
  //char *fileptr = new char [infile.length()+1];
  //std::strcpy(fileptr, infile.c_str());

  // char *locatedFile =  prefixFilename(fileptr, flukarun);
  // std::cout << "prefixed file is " << locatedFile << std::endl;
//  return checkFile(locatedFile);
//}


//---------------------------------------------------------------------------//
// checkFile(..)
//---------------------------------------------------------------------------//
/// Call inline function to check basic i/o
bool checkFile(char *fileptr)
{
     std::cerr << __FILE__ << ", " << __func__ << ":" << __LINE__ << "_______________" << std::endl;
     std::cerr << "Checking file " << fileptr << std::endl;
     bool readable = isFileReadable(fileptr);
     std::cerr << "Is " << fileptr << " readable?  " << readable << std::endl;
     bool nonEmpty = isFileNonEmpty(fileptr);
     std::cerr << "Is " << fileptr << " nonempty?  " << nonEmpty << std::endl;
     bool both = readable && nonEmpty;
     std::cerr << "Is " << fileptr << " both?  " << both << std::endl;

     return both;

//   return isFileReadable(fileptr) && isFileNonEmpty(fileptr); 
}

bool checkFile(std::string filename)
{
  char *fileptr = new char[filename.length() + 1];
  std::strcpy(fileptr, filename.c_str());
 // return file_can_be_read(fileptr); 
  bool success = std::ifstream(fileptr);
  if (success)
  {
    std::cerr << "Input file: " << filename << " can be read." << std::endl;
  }
  return success;
}


//---------------------------------------------------------------------------//
// checkFile(..)
//---------------------------------------------------------------------------//
// Test the non-prefixed file because the testing is done prior to calling fluka
bool isFileReadable(char *fileptr)
{
  
  bool success =   std::ifstream(fileptr); 
  if (success)
  {
    std::cerr << "Input file: " << fileptr << " can be read." << std::endl;
  }
  return success;
}

//---------------------------------------------------------------------------//
// isFileNonEmpty(..)
//---------------------------------------------------------------------------//
// Test the prefixed file because the testing is done prior to calling fluka
bool isFileNonEmpty(char *fileptr)
{
	// std::string s = fileptr;
	std::ifstream file(fileptr, std::ios::binary);
        
	file.seekg(0, std::ios::end);
	if (file.tellg() == 0)
	{
		std::cerr << "File is empty" << std::endl;
                return false;
	}
        else
        {
               return true;
        }
}
