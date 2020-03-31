// Source.cpp
// Central Source Class

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

pyne::PointSource::PointSource(double x, double y, double z, double u, double v,
                               double k, double E, std::string p, double weight)
    : x(x), y(y), z(z), i(i), j(j), k(k), E(E), p(p), weight(weight) {}

// Destructor
pyne::PointSource::~PointSource() {}

std::string pyne::PointSource::mcnp() {
  std::stringstream output;  // output stream
  std::string indent = "     ";

  output << "SDEF ";

  if (x != 0 && y != 0 && z != 0) {
    output << "POS=" << x << " " << y << " " << z << std::endl;
  }

  if (i != 0 && j != 0 && k != 0) {
    output << indent << "VEC=" << i << " " << j << " " << k << "DIR=1"
           << std::endl;
  }

  if (E != 14) {
    output << std::endl << indent << "ERG=" << E;
  }

  if (weight != 1) {
    output << std::endl << indent << "WGT=" << weight;
  }

  if (p != "n") {
    output << std::endl << indent << "PAR=" << pyne::particle::mcnp(p);
  }

  return output.str();
}

// Produces valid fluka source
std::string pyne::PointSource::fluka() {
  std::stringstream output;  // output stream

  return output.str();
}
