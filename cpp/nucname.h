// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// name is for letters  as well (U-235).
// MCNP is for numerals without the meta-stable flag (92235), as used in MCNP.

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


  /// \name Current Form Functions
  /// \{
  /// These functions compute the current form of the nuclide.
  /// \param nuc a nuclide
  /// \return a (string) flag of the nuclides current form.
  std::string current_form(std::string nuc);
  std::string current_form(int nuc);
  /// \}

  /// \name ZZAAAM Form Functions
  /// \{
  /// The ZZAAAM nuclide naming convention is the canonical form for representing
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


  /************************/
  /*** name functions ***/
  /************************/
  std::string name(int);
  std::string name(char *);
  std::string name(std::string);


  /**********************/
  /*** mcnp functions ***/
  /**********************/
  int mcnp(int);
  int mcnp(char *);
  int mcnp(std::string);


  /*************************/
  /*** serpent functions ***/
  /*************************/
  std::string serpent(int);
  std::string serpent(char *);
  std::string serpent(std::string);


  /**********************/
  /*** nist functions ***/
  /**********************/
  std::string nist(int);
  std::string nist(char *);
  std::string nist(std::string);


  /************************/
  /*** cinder functions ***/
  /************************/
  int cinder(int);
  int cinder(char *);
  int cinder(std::string);

};
};

#endif
