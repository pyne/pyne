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
  virtual std::string mcnp();

  // fluka tally
  virtual std::string fluka();
}

class SourcePoint : public Source {
 public:
  /// Point Source Constructors
  SourcePoint(double x = 0, double y = 0, double z = 0, double u = 0,
              double v = 0, double w = 0, double E = 14, std::string p = "n",
              double weight = 1);  /// default constructor

  ~SourcePoint();  /// default destructor

  // mcnp tally
  std::string mcnp();

  // fluka tally
  std::string fluka();

 private:
  double x;
  double y;
  double z;

  double u;
  double v;
  double w;

  double E;
  double weight;
  std::string p;
};

/// Converts a Source to a string stream representation.
std::ostream& operator<<(std::ostream& os, Source source);

// End pyne namespace
}  // namespace pyne

#endif
