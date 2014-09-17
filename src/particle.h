/// \brief Converts betweeen naming/numbering conventions for particle types

// defines the primary particle types that are allowed by most monte carlo codes
// some monte carlo codes allow us to score so called "heavy ions", in fact we
// define heavy ions to be particles with more than one neutron or proton

#include <string>
#include <map>
#include <set>

#ifndef PYNE_IS_AMALGAMATED
  #include "nucname.h"
#endif

/// Number of pure particle types currently supported
#define NUM_PARTICLES 32

namespace pyne 
{
namespace particle
{
  extern int _pdcids[NUM_PARTICLES]; ///
  /// set of Particle Data Centre integer id numbers
  extern std::string _docs[NUM_PARTICLES];
  /// set of doc strings that describe the particle types
  extern std::string _names[NUM_PARTICLES];
  /// set of name strings that are the particle types
  extern std::set<std::string> names;
  /// set of valid names
  extern std::set<int> pdc_nums;
  /// set of valid pdc numbers
  extern std::map<std::string,int> name_id;
  /// map of name to pdc number
  extern std::map<int,std::string> id_name;
  /// map of pdc number to name
  extern std::map<std::string,std::string> docs;
  /// map of name to doc string
  extern std::map<std::string,int> altnames;
  /// map of alternative name to pdc number

  /// \name is_hydrogen functions
  /// \{
  /// Returns whether or not the given particle is hydrogen or not, for
  /// example, Protons (Hydrogen) are both valid nucids and fundamental 
  /// pdc types, all the following identify as hydrogen, Proton, Hydrogen,
  /// Protium, "H1", "1H", 100001000, PDC(2212)
  /// \param n Integer PDC number or nucid
  /// \param s String valid particle name, altname or nucid
  bool _is_hydrogen(int n);
  bool _is_hydrogen(char *s);
  bool _is_hydrogen(std::string s);
  /// \}

  /// \name is_heavy_ion functions
  /// \{
  /// Returns whether or not the given particle is a heavy ion
  /// or not. Heavy ions are not covered by the PDC scheme, therefore
  /// the pyne::nucname class is used.
  /// \param n Integer PDC number or nucid
  /// \param s String valid particle name, altname or nucid
  bool is_heavy_ion(int s);
  bool is_heavy_ion(char *s);
  bool is_heavy_ion(std::string s);
  ///\}

  /// \name is_valid functions
  /// \{
  /// Returns whether or not the the given particle is a valid particle
  /// in the PyNE particle class. All PDC numbers, names, altnames nucnames
  /// are valid particle types
  /// \param n Integer PDC number or nucid
  /// \param s String valid particle name, altname or nucid
  bool is_valid(int n);
  bool is_valid(char *s);
  bool is_valid(std::string s);
  ///\}

  /// \name pdc_number functions
  /// \{
  /// Returns the PDC number of the particle given, if a valid pdc particle, 
  /// will return the number, for heavy ions will return 0.
  /// \param n Integer PDC number or nucid
  /// \param s String valid particle name, altname or nucid
  int pdc_number(int n);
  int pdc_number(char *s);
  int pdc_number(std::string s);
  ///\}

  /// \name name functions
  /// \{
  /// Returns the pyne::particle name of the particle given, 
  /// if a valid pdc particle number, name or nucname. Raises
  /// exception if not a valid name
  /// \param n Integer PDC number or nucid
  /// \param s String valid particle name, altname or nucid
  std::string name(int s);
  std::string name(char *s);
  std::string name(std::string s);
  ///\}

  /// \name describe functions
  /// \{
  /// Returns a long string that describes the particle, if given
  /// a valid particle name, otherwise raises exception
  /// \param n Integer PDC number or nucid
  /// \param s String valid particle name, altname or nucid
  std::string describe(int s);
  std::string describe(char *s);
  std::string describe(std::string s);

  /// A helper function to set the contents of the variables in this library.
  void * _fill_maps();  
  extern void * _;  ///< A dummy variable used when calling #_fill_maps().

  /// Custom excpeption for failed particle types
  class NotAParticle : public std::exception
  {
  public:
    /// Default constructor
    NotAParticle () {};

    /// Default destructor
    ~NotAParticle () throw () {};

    /// Constructor for raising the exception
    /// Spits out the particle name as input
    NotAParticle(std::string particle_name)
    {
       part_name = particle_name;
    };

    /// raises error message
    virtual const char* what() const throw()
    {
      std::string pname ("Not a valid particle name ");
      if(!part_name.empty())
	pname += part_name;
      return (const char *) pname.c_str();
    };

    private:
       std::string part_name;  /// the particle name


  };
};
};
  
