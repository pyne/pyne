#include <fstream>

inline bool file_exists_and_can_be_read_from( const char* path2file )
{ 
    return std::ifstream(path2file); 
}

// Handle logic and errors in the input
char* checkInput(std::string infile, bool& flukarun);

// Check for existence and readabilieyt of a file, handling failure
bool checkFile(char* fileptr);
bool isFileReadable(char *fileptr);
bool isFileNonEmpty(char *fileptr);

// When passing a filename to be used as fluka input, it must be prefixed
// with "../" because it will be referenced from a subdirectory.
char *prefixFilename(char *cfile, bool running_with_fluka);
