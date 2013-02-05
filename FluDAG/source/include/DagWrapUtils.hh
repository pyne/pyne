
// after flugg Wrappers/include/WrapUtils.hh

///////////////////////////////////////////////////////////////////
//
// Some utility methods used for several wrapped functions
//
///////////////////////////////////////////////////////////////////

#ifndef DagWrapUtils_hh
#define DagWrapUtils_hh 1

#include <iomanip>

//Forward declarations

// Common print utilities used in FGeometryInit.cc
inline std::ostream& setw10(std::ostream& os) { return os << std::setw(10);}
inline std::ostream& setscientific(std::ostream& os) { return os << std::setiosflags(std::ios::scientific);}
inline std::ostream& setfixed(std::ostream& os) { return os << std::setiosflags(std::ios::fixed);}
std::ostream& PrintHeader(std::ostream& os, const char* title);

#endif
