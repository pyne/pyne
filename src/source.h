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


namespace pyne
{
  // simplifying implementation...
  using std::vector;
  // default empty parameters for constructors
  const vector<double> null_v_dbl;
  const vector<int> null_v_int;

  class Source 
  {
  public:

    /// Source Constructors
    Source (); /// empty constructor


    ~Source();  /// default destructor
 

    // mcnp tally
    std::string mcnp();
    
    // fluka tally
    std::string fluka(std::string unit_number = "-21");

  };

  /// Converts a Source to a string stream representation.
  std::ostream& operator<< (std::ostream& os, Source source);



// End pyne namespace
}

#endif
