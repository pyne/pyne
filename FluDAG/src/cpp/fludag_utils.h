//----------------------------------*-C++, Fortran-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/include/fludag_utils.h
 * \author Julie Zachman 
 * \date   Apr 5 2013 
 * \brief  Functions called from Input.cpp and mainFluDAG.cpp
 * \note   Unittested
 */
//---------------------------------------------------------------------------//
#include <fstream>

// If you can create an input filestream with it, it's there.
inline bool file_can_be_read( const char* path2file )
{ 
    return std::ifstream(path2file); 
}

// Perform the task of argument and error checking
void checkArgs(std::string infile, int numargs, char* argv[], bool& flukarun);
// Use flag to decide the relative location of input files
std::string prefixFilename(std::string infile, bool running_with_fluka);

// Handle logic and errors in the input
// Check for existence and readability of a file, handling failure
// Pass in the string and create a const char * copy from it.
bool checkFile(std::string filename);
// Call the inline above with correct arg
bool isFileReadable(std::string filename);
// Make sure the file is not of length zero
bool isFileNonEmpty(std::string filename);

