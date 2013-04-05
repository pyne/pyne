#include "fludag_utils.h"
#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string.h>

#define DEBUG 1

// Precondition: argument number has already been checked
// Handle error conditions gracefully
char* checkInput(std::string infile, bool& flukarun)
{
  char *fileptr = new char [infile.length()+1];
  std::strcpy(fileptr, infile.c_str());

  if (checkFile(fileptr))
  {
    return prefixFilename(fileptr, flukarun);
  }
  else
  {
    return NULL;
  }
}


bool checkFile(char *fileptr)
{
     return isFileReadable(fileptr) && isFileNonEmpty(fileptr); 
}

// Test the non-prefixed file because the testing is done prior to calling fluka
bool isFileReadable(char *fileptr)
{
  if ( file_exists_and_can_be_read_from(fileptr) )
  {
     return true;
  }
  {
     std::cerr << "Input file: " << fileptr << " cannot be read." << std::endl;
     return false;
  }
}

bool isFileNonEmpty(char *fileptr)
{
	// std::string s = fileptr;
	std::ifstream file(fileptr, std::ios::binary);
        
	file.seekg(0, std::ios::end);
	if (file.tellg() == 0)
	{
		std::cout << "File is empty" << std::endl;
                return false;
	}
        else
        {
               return true;
        }
}

// If running_with_fluka, prefix "../" to the filename, as fluka runs from a subdirectory
char *prefixFilename(char *cfile, bool running_with_fluka)
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
  	std::cout << "\nmy file is " << prefixedFilename << "\n" << std::endl;
  }
  char *myfile;
  myfile = &prefixedFilename[0];
  return myfile;
}  


