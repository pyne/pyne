// Source.cpp
// Central Source Class
// -- Andrew Davis

#include <iomanip>
#include <string>
#include <vector>

#ifndef PYNE_IS_AMALGAMATED
  #include "particle.h"
  #include "source.h"
#endif


/***************************/
/*** Protected Functions ***/
/***************************/

// there are no protected functions currently
// fool.

/************************/
/*** Public Functions ***/
/************************/

/*--- Constructors ---*/

pyne::Source::Source() {
}

// Destructor
pyne::Source::~Source() {}

std::ostream& operator<<(std::ostream& os, pyne::Source tal) {
  //print the Source to ostream
  os << "\t---------\n";
  return os;
}

// Sets string to valid mcnp formatted source
// Takes mcnp version as arg, like 5 or 6
std::string pyne::Source::mcnp() {
  std::stringstream output;  // output stream
  
  return output.str();
}

// Produces valid fluka source
std::string pyne::Source::fluka() {
  std::stringstream output; // output stream
  
  return output.str();
}
