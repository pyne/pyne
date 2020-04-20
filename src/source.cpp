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

pyne::PointSource::PointSource(double _x, double _y, double _z, double _u, double _v,
                               double _w, double _E, std::string _particle, double _weight)
    : x(_x), y(_y), z(_z), u(_u), v(_v), w(_w), E(_E), particle(_particle), weight(_weight) {}

// Destructor
pyne::PointSource::~PointSource() {}

std::string pyne::PointSource::mcnp(int version) const {
  std::stringstream output;  // output stream

  std::stringstream newline;
  newline << std::endl << "     ";

  output << "SDEF ";
  output << "POS=" << x << " " << y << " " << z;
  if (u != 0 || v != 0 || w != 0) {
    output << newline.str();
    output << "VEC=" << u << " " << v << " " << w << " DIR=1";
  }
  output << newline.str() << "ERG=" << E;
  output << newline.str() << "WGT=" << weight;
  output << newline.str() << "PAR=";
  if (version == 5)
    output << pyne::particle::mcnp(particle);
  else if (version == 6)
    output << pyne::particle::mcnp6(particle);
  else {
    throw std::runtime_error("MCNP version is not reconized!");
  }

  return output.str();
}

// Produces valid fluka source
std::string pyne::PointSource::fluka() {
  std::stringstream output;  // output stream

  return output.str();
}
