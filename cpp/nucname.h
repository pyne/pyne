/// \file nucname.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Converts between naming conventions for nuclides.

#ifndef _PYNE_NUCNAME_
#define _PYNE_NUCNAME_
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>

#include "pyne.h"

namespace pyne
{
namespace nucname 
{
  typedef std::string name_t; ///< name type
  typedef int zz_t;           ///< Z number type

  typedef std::map<name_t, zz_t> name_zz_t; ///< name and Z num map type
  typedef name_zz_t::iterator name_zz_iter; ///< name and Z num iter type
  name_zz_t get_name_zz();  ///< Creates standard name to Z number mapping.
  extern name_zz_t name_zz; ///< name to Z num map

  typedef std::map<zz_t, name_t> zzname_t;  ///< Z num to name map type
  typedef zzname_t::iterator zzname_iter;   ///< Z num to name iter type
  zzname_t get_zz_name();   ///< Creates standard Z number to name mapping.
  extern zzname_t zz_name;  ///< Z num to name map

  /******************************************/
  /*** Define useful elemental group sets ***/
  /******************************************/

  /// name grouping type (for testing containment)
  typedef std::set<name_t> name_group;  
  typedef name_group::iterator name_group_iter; ///< name grouping iter type

  /// Z number grouping type (for testing containment)
  typedef std::set<zz_t> zz_group;
  typedef zz_group::iterator zz_group_iter; ///< Z number grouping iter

  /// Converts a name group to a Z number group.
  /// \param eg a grouping of nuclides by name
  /// \return a Z numbered group
  zz_group name_to_zz_group (name_group eg);

  extern name_t LAN_array[15];  ///< array of lanthanide names
  extern name_group LAN;        ///< lanthanide name group
  extern zz_group lan;          ///< lanthanide Z number group

  extern name_t ACT_array[15];  ///< array of actinide names
  extern name_group ACT;        ///< actinide name group
  extern zz_group act;          ///< actinide Z number group

  extern name_t TRU_array[22];  ///< array of transuranic names
  extern name_group TRU;        ///< transuranic name group
  extern zz_group tru;          ///< transuranic Z number group

  extern name_t MA_array[10];   ///< array of minor actinide names
  extern name_group MA;         ///< minor actinide name group
  extern zz_group ma;           ///< minor actinide Z number group

  extern name_t FP_array[88];   ///< array of fission product names
  extern name_group FP;         ///< fission product name group
  extern zz_group fp;           ///< fission product Z number group


  /******************/
  /*** Exceptions ***/
  /******************/

  /// Custom expection for declaring that a value does not follow a recognizable 
  /// nuclide naming convention.
  class NotANuclide : public std::exception
  {
  public:
    /// default constructor
    NotANuclide () {}; 

    /// default destructor
    ~NotANuclide () throw () {};

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    NotANuclide(std::string wasptr, std::string nowptr)
    {
       nucwas = wasptr;
       nucnow = nowptr;
    };

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    NotANuclide(std::string wasptr, int nowptr)
    {
      nucwas = wasptr;
      nucnow = pyne::to_str(nowptr);
    };

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    NotANuclide(int wasptr, std::string nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = nowptr;
    };

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    NotANuclide(int wasptr, int nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = pyne::to_str(nowptr);
    };

    /// Generates an informational message for the exception
    /// \return The error string
    virtual const char* what() const throw()
    {
      std::string NaNEstr ("Not a Nuclide! ");
      if (!nucwas.empty())
        NaNEstr += nucwas;

      if (!nucnow.empty())
      {
        NaNEstr += " --> "; 
        NaNEstr += nucnow;
      }
      return (const char *) NaNEstr.c_str();
    };

  private:
    std::string nucwas; ///< previous nuclide state
    std::string nucnow; ///< current nuclide state  
  };

  /// Custom expection for declaring that a value represents one or more nuclides 
  /// in one or more namig conventions
  class IndeterminateNuclideForm : public std::exception
  {
  public:
    /// default constructor
    IndeterminateNuclideForm () {};

    /// default destuctor
    ~IndeterminateNuclideForm () throw () {};

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    IndeterminateNuclideForm(std::string wasptr, std::string nowptr)
    {
       nucwas = wasptr;
       nucnow = nowptr;
    };

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    IndeterminateNuclideForm(std::string wasptr, int nowptr)
    {
      nucwas = wasptr;
      nucnow = pyne::to_str(nowptr);
    };

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    IndeterminateNuclideForm(int wasptr, std::string nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = nowptr;
    };

    /// Constructor given previous and current state of nulide name
    /// \param wasptr Previous state, typically user input.
    /// \param nowptr Current state, as far as PyNE could get.
    IndeterminateNuclideForm(int wasptr, int nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = pyne::to_str(nowptr);
    };

    /// Generates an informational message for the exception
    /// \return The error string
    virtual const char* what() const throw()
    {
      std::string INFEstr ("Indeterminate nuclide form: ");
      if (!nucwas.empty())
        INFEstr += nucwas;

      if (!nucnow.empty())
      {
        INFEstr += " --> "; 
        INFEstr += nucnow;
      }
      return (const char *) INFEstr.c_str();
    }

  private:
    std::string nucwas; ///< previous nuclide state
    std::string nucnow; ///< current nuclide state
  };

  /// \name isnuclide functions
  /// \{
  /// These functions test if an input \a nuc is a valid nuclide.
  /// \param nuc a possible nuclide
  /// \return a bool
  bool isnuclide(std::string nuc);
  bool isnuclide(char * nuc);
  bool isnuclide(int nuc);
  /// \}

  /// \name Identifier Form Functions
  /// \{
  /// The 'id' nuclide naming convention is the canonical form for representing
  /// nuclides in PyNE. This is termed a ZAS, or ZZZAAASSSS, representation because 
  /// It stores 3 Z-number digits, 3 A-number digits, followed by 4 S-number digits
  /// which the nucleus excitation state. 
  /// \param nuc a nuclide
  /// \return nucid 32-bit integer identifier
  int id(int nuc);
  int id(char * nuc);
  int id(std::string nuc);
  /// \}

  /// \name Name Form Functions
  /// \{
  /// The 'name' nuclide naming convention is the more common, human readable 
  /// notation. The chemical symbol (one or two characters long) is first, followed 
  /// by the nucleon number. Lastly if the nuclide is metastable, the letter M is 
  /// concatenated to the end. For example, ‘H-1’ and ‘Am242M’ are both valid. 
  /// Note that nucname will always return name form with the dash removed and all 
  /// letters uppercase.
  /// \param nuc a nuclide
  /// \return a string nuclide identifier.
  std::string name(int nuc);
  std::string name(char * nuc);
  std::string name(std::string nuc);
  /// \}

  /// \name Z-Number Functions
  /// \{
  /// The Z-number, or charge number, represents the number of protons in a
  /// nuclide.  This function returns that number.
  /// \param nuc a nuclide
  /// \return an integer Z-number.
  int znum(int nuc);
  int znum(char * nuc);
  int znum(std::string nuc);
  /// \}

  /// \name A-Number Functions
  /// \{
  /// The A-number, or nucleon number, represents the number of protons and 
  /// neutrons in a nuclide.  This function returns that number.
  /// \param nuc a nuclide
  /// \return an integer A-number.
  int anum(int nuc);
  int anum(char * nuc);
  int anum(std::string nuc);
  /// \}

  /// \name S-Number Functions
  /// \{
  /// The S-number, or excitation state number, represents the excitation
  /// level of a nuclide.  Normally, this is zero.  This function returns 
  /// that number.
  /// \param nuc a nuclide
  /// \return an integer A-number.
  int snum(int nuc);
  int snum(char * nuc);
  int snum(std::string nuc);
  /// \}

  /// \name ZZAAAM Form Functions
  /// \{
  /// The ZZAAAM nuclide naming convention is the former canonical form for 
  /// nuclides in PyNE. This places the charge of the nucleus out front, then has 
  /// three digits for the atomic mass number, and ends with a metastable flag 
  /// (0 = ground, 1 = first excited state, 2 = second excited state, etc). 
  /// Uranium-235 here would be expressed as ‘922350’.
  /// \param nuc a nuclide
  /// \return an integer nuclide identifier.
  int zzaaam(int nuc);
  int zzaaam(char * nuc);
  int zzaaam(std::string nuc);
  /// \}

  /// \name ZZAAAM Form to Identifier Form Functions
  /// \{
  /// This converts from the ZZAAAM nuclide naming convention 
  /// to the id canonical form  for nuclides in PyNE. 
  /// \param nuc a nuclide in ZZAAAM form.
  /// \return an integer id nuclide identifier.
  int zzaaam_to_id(int nuc);
  int zzaaam_to_id(char * nuc);
  int zzaaam_to_id(std::string nuc);
  /// \}

  /// \name MCNP Form Functions
  /// \{
  /// This is the naming convention used by the MCNP suite of codes.
  /// The MCNP format for entering nuclides is unfortunately non-standard. 
  /// In most ways it is similar to zzaaam form, except that it lacks the metastable 
  /// flag. For information on how metastable isotopes are named, please consult the 
  /// MCNP documentation for more information.
  /// \param nuc a nuclide
  /// \return a string nuclide identifier.
  int mcnp(int nuc);
  int mcnp(char * nuc);
  int mcnp(std::string nuc);
  /// \}

  /// \name MCNP Form to Identifier Form Functions
  /// \{
  /// This converts from the MCNP nuclide naming convention 
  /// to the id canonical form  for nuclides in PyNE. 
  /// \param nuc a nuclide in MCNP form.
  /// \return an integer id nuclide identifier.
  int mcnp_to_id(int nuc);
  int mcnp_to_id(char * nuc);
  int mcnp_to_id(std::string nuc);
  /// \}

  /// \name Serpent Form Functions
  /// \{
  /// This is the string-based naming convention used by the Serpent suite of codes.
  /// The serpent naming convention is similar to name form. However, only the first 
  /// letter in the chemical symbol is uppercase, the dash is always present, and the 
  /// the meta-stable flag is lowercase. For instance, ‘Am-242m’ is the valid serpent 
  /// notation for this nuclide.
  /// \param nuc a nuclide
  /// \return a string nuclide identifier.
  std::string serpent(int nuc);
  std::string serpent(char * nuc);
  std::string serpent(std::string nuc);
  /// \}

  /// \name Serpent Form to Identifier Form Functions
  /// \{
  /// This converts from the Serpent nuclide naming convention 
  /// to the id canonical form  for nuclides in PyNE. 
  /// \param nuc a nuclide in Serpent form.
  /// \return an integer id nuclide identifier.
  //int serpent_to_id(int nuc);  Should be ZAID
  int serpent_to_id(char * nuc);
  int serpent_to_id(std::string nuc);
  /// \}

  /// \name NIST Form Functions
  /// \{
  /// This is the string-based naming convention used by NIST.
  /// The NIST naming convention is also similar to the Serpent form. However, this 
  /// convention contains no metastable information. Moreover, the A-number comes 
  /// before the element symbol. For example, ‘242Am’ is the valid NIST notation.
  /// \param nuc a nuclide
  /// \return a string nuclide identifier.
  std::string nist(int nuc);
  std::string nist(char * nuc);
  std::string nist(std::string nuc);
  /// \}

  /// \name NIST Form to Identifier Form Functions
  /// \{
  /// This converts from the NIST nuclide naming convention 
  /// to the id canonical form  for nuclides in PyNE. 
  /// \param nuc a nuclide in NIST form.
  /// \return an integer id nuclide identifier.
  //int serpent_to_id(int nuc);  NON-EXISTANT
  int nist_to_id(char * nuc);
  int nist_to_id(std::string nuc);
  /// \}

  /// \name CINDER Form Functions
  /// \{
  /// This is the naming convention used by the CINDER burnup library.
  /// The CINDER format is similar to zzaaam form except that the placement of the 
  /// Z- and A-numbers are swapped. Therefore, this format is effectively aaazzzm. 
  /// For example, ‘2420951’ is the valid cinder notation for ‘AM242M’.
  /// \param nuc a nuclide
  /// \return a string nuclide identifier.
  int cinder(int nuc);
  int cinder(char * nuc);
  int cinder(std::string nuc);
  /// \}

  /// \name Cinder Form to Identifier Form Functions
  /// \{
  /// This converts from the Cinder nuclide naming convention 
  /// to the id canonical form  for nuclides in PyNE. 
  /// \param nuc a nuclide in Cinder form.
  /// \return an integer id nuclide identifier.
  int cinder_to_id(int nuc);
  int cinder_to_id(char * nuc);
  int cinder_to_id(std::string nuc);
  /// \}

  /// \name ALARA Form Functions
  /// \{
  /// This is the format used in the ALARA activation code elements library.
  /// For elements, the form is "ll" where ll is the atomic symbol. For isotopes
  /// the form is "ll:AAA". No metastable isotope flag is used.
  /// \param nuc a nuclide
  /// \return a string nuclide identifier.
  std::string alara(int nuc);
  std::string alara(char * nuc);
  std::string alara(std::string nuc);
  /// \}

  /// \name ALARA Form to Identifier Form Functions
  /// \{
  /// This converts from the ALARA nuclide naming convention 
  /// to the id canonical form  for nuclides in PyNE. 
  /// \param nuc a nuclide in ALARA form.
  /// \return an integer id nuclide identifier.
  //int alara_to_id(int nuc); NOT POSSIBLE
  int alara_to_id(char * nuc);
  int alara_to_id(std::string nuc);
  /// \}

  /// \name SZA Form Functions
  /// \{
  /// This is the new format for ACE data tables in the form SSSZZZAAA.
  /// The first three digits represent the excited state (000 = ground,
  /// 001 = first excited state, 002 = second excited state, etc).
  /// The second three digits are the atomic number and the last three
  /// digits are the atomic mass. Prepending zeros can be omitted, making
  /// the SZA form equal to the MCNP form for non-excited nuclides.
  /// \param nuc a nuclide
  /// \return a string nuclide identifier.
  int sza(int nuc);
  int sza(char * nuc);
  int sza(std::string nuc);
  /// \}

  /// \name SZA Form to Identifier Form Functions
  /// \{
  /// This converts from the SZA nuclide naming convention 
  /// to the id canonical form  for nuclides in PyNE. 
  /// \param nuc a nuclide in SZA form.
  /// \return an integer id nuclide identifier.
  int sza_to_id(int nuc);
  int sza_to_id(char * nuc);
  int sza_to_id(std::string nuc);
  /// \}

};
};

#endif
