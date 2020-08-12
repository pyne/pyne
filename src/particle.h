/// \brief Converts betweeen naming/numbering conventions for particle types

// defines the primary particle types that are allowed by most monte carlo codes
// some monte carlo codes allow us to score so called "heavy ions", in fact we
// define heavy ions to be particles with more than one neutron or proton

#ifndef PYNE_UWPOP4EE6BEB5CZK4BQQHYTCEI
#define PYNE_UWPOP4EE6BEB5CZK4BQQHYTCEI

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
//! Converts betweeen naming/numbering conventions for particle types
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
  extern std::map<std::string,std::string> part_to_mcnp;
  /// map of name to mcnp string
  extern std::map<std::string,std::string> part_to_mcnp6;
  /// map of name to mcnp6 string
  extern std::map<std::string,std::string> part_to_fluka;
  /// map of name to fluka string
  extern std::map<std::string,std::string> part_to_geant4;
  /// map of name to geant4 string


  /// \name is_hydrogen functions
  /// \{
  /// Returns whether or not the given particle is hydrogen or not, for
  /// example, Protons (Hydrogen) are both valid nucids and fundamental
  /// pdc types, all the following identify as hydrogen, Proton, Hydrogen,
  /// Protium, "H1", "1H", 100001000, PDC(2212)
  /// \param n Integer PDC number or nucid
  /// \param s String valid particle name, altname or nucid
  bool is_hydrogen(int n);
  bool is_hydrogen(char *s);
  bool is_hydrogen(std::string s);
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
  int id(int n);
  int id(char *s);
  int id(std::string s);
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

  /// \name mcnp functions
  /// \{
  /// Returns the mcnp string of a valid pyne::particle name
  /// \param s int, char*, String valid particle name, altname or nucid
  std::string mcnp(int s);
  std::string mcnp(char *s);
  std::string mcnp(std::string s);
  ///\}

  /// \name mcnp6 functions
  /// \{
  /// Returns the mcnp6 string of a valid pyne::particle name
  /// \param s int, char*, String valid particle name, altname or nucid
  std::string mcnp6(int s);
  std::string mcnp6(char *s);
  std::string mcnp6(std::string s);
  ///\}

  /// \name fluka functions
  /// \{
  /// Returns the Fluka string of a valid pyne::particle name, or heavy ion
  /// \param s int, char*, String valid particle name, altname or nucid
  std::string fluka(int s);
  std::string fluka(char *s);
  std::string fluka(std::string s);
  ///\}

  /// \name geant4 functions
  /// \{
  /// Returns the Geant4 string of a valid pyne::particle name, or heavy ion
  /// \param s int, char*, String valid particle name, altname or nucid
  std::string geant4(int s);
  std::string geant4(char *s);
  std::string geant4(std::string s);
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
  extern void * filler;  ///< A dummy variable used when calling #_fill_maps().


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
    }

    /// raises error message
    virtual const char* what() const throw()
    {
      std::string pname ("Not a valid particle name ");
      if(!part_name.empty())
        pname += part_name;
      const char *pname_rtn = pname.c_str();
      return pname_rtn;
    }

    private:
       std::string part_name;  /// the particle name


  };
}
}

#endif
