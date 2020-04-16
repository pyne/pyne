/// \brief The tally class and helper functions.
///
/// The tally class is in essesence a structure containing attributes
/// related to tallies.

#ifndef PYNE_IQ4M73STINHJDPRV6KWUZZXOYE
#define PYNE_IQ4M73STINHJDPRV6KWUZZXOYE

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#ifndef PYNE_IS_AMALGAMATED
#include "h5wrap.h"
#include "utils.h"
#endif

namespace pyne {

class Source {
 public:
  /// Source Constructors
  Source();  /// empty constructor

  ~Source();  /// default destructor

  // mcnp tally
  virtual std::string mcnp(int version=5) const {return "";}

  // fluka tally
  virtual std::string fluka() {}
};

class PointSource : public Source {
 public:
  /// Point Source Constructors
  PointSource(double _x = 0, double _y = 0, double _z = 0, double _u = 0,
              double _v = 0, double _w = 0, double _E = 14, std::string _particle = "Neutron",
              double _weight = 1);  /// default constructor

  ~PointSource();  /// default destructor

  // mcnp tally
  virtual std::string mcnp(int version=5) const;

  // fluka tally
  std::string fluka();



  double x;
  double y;
  double z;

  double u;
  double v;
  double w;

  double E;
  double weight;
  std::string particle;
};

/// Converts a Source to a string stream representation.
std::ostream& operator<<(std::ostream& os, Source source);

// End pyne namespace
}  // namespace pyne

#endif
