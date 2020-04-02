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

pyne::PointSource::PointSource(double x, double y, double z, double i, double j,
                               double k, double E, std::string p, double weight)
    : x(x), y(y), z(z), i(i), j(j), k(k), E(E), particle(p), weight(weight) {}

// Destructor
pyne::PointSource::~PointSource() {}

std::string pyne::PointSource::mcnp(int version) {
  std::stringstream output;  // output stream
  
  std::stringstream newline;
  newline << std::endl << "     ";
  
  
  output << "SDEF ";
  if (x != 0 || y != 0 || z != 0) {
    output << "POS=" << x << " " << y << " " << z;
  }
  if (i != 0 || j != 0 || k != 0) {
    output << newline.str();
    output << "VEC=" << i << " " << j << " " << k << " DIR=1";
  }
  if (E != 14) {
    output << newline.str() << "ERG=" << E;
  }
  if (weight != 1) {
    output << newline.str() << "WGT=" << weight;
  }
  if (particle != "n") {
    output << newline.str() << "PAR=";
    if (version == 5)
      output << pyne::particle::mcnp(particle);
    else if (version == 6)
      output << pyne::particle::mcnp6(particle);
    else {
      throw std::runtime_error("MCNP version is not reconized!");
    }
  }

  return output.str();
}

// Produces valid fluka source
std::string pyne::PointSource::fluka() {
  std::stringstream output;  // output stream

  return output.str();
}
