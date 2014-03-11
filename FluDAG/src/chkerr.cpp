#include "chkerr.hpp"

using namespace moab;

void chkerr( Interface* mbi, ErrorCode code, int line, const char* file )
{
  if( code != MB_SUCCESS )
  {
    std::cerr << "Error: " << mbi->get_error_string(code) << " (" << code << ")" << std::endl;
    std::string message;
    if (MB_SUCCESS == mbi->get_last_error(message) && !message.empty())
    {
      std::cerr << "Error message: " << message << std::endl;
    }
    std::string fname = file;
    size_t slash = fname.find_last_of('/');
    if( slash != fname.npos ){ fname = fname.substr(slash+1); }
    std::cerr << "At " << fname << " line: " << line  << std::endl;
    std::exit( EXIT_FAILURE );
  }
}

void chkerr( Interface& mbi, ErrorCode code, int line, const char* file )
{
  chkerr( &mbi, code, line, file );
}
 
