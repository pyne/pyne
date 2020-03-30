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

pyne::Source::Source() {}

// Destructor
pyne::Source::~Source() {}

std::ostream& operator<<(std::ostream& os, pyne::Source tal) {
  // print the Source to ostream
  os << "\t---------\n";
  return os;
}

pyne::SourcePoint::SourcePoint(double x, double y, double z, double u, double v,
                               double w, double E, std::string p, double weight)
    : x(x), y(y), z(z), u(u), v(v), w(w), E(E), p(p), weight(weight) {}

    
// Destructor
pyne::SourcePoint::~SourcePoint() {}

std::string pyne::SourcePoint::mcnp() {
  std::stringstream output;  // output stream
  std::string separator = "     "

                          output
                          << "SDEF "

                          if (x != 0 && y != 0 && z != 0) {
    output << "POS=" << x << " " << y << " " << z << std::endl;
  }

  if (E != 14) {
    output << std::endl << separator << "ERG=" << E;
  }

  if (weight != 1) {
    output << std::endl << separator << "WGT=" << weight;
  }

  if (p != "n") {
    output << std::endl << separator << "PAR=" << pyne::particle::mcnp(p);
  }

  return output.str();
}

// Produces valid fluka source
std::string pyne::SourcePoint::fluka() {
  std::stringstream output;  // output stream

  return output.str();
}
