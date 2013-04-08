#include <fstream>

inline bool file_can_be_read( const char* path2file )
{ 
    return std::ifstream(path2file); 
}

// Handle logic and errors in the input
bool checkInput(char *fileptr);

// Check for existence and readability of a file, handling failure
// bool checkFile(char* fileptr);
// Pass in the string and create a const char * copy from it.
bool checkFile(std::string filename);
bool isFileReadable(char *fileptr);
bool isFileNonEmpty(char *fileptr);

// When passing a filename to be used as fluka input, it must be prefixed
// with "../" because it will be referenced from a subdirectory.
// char *prefixFilename(char *cfile, bool running_with_fluka);
