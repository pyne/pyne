/// \file rxname.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Converts between naming conventions for reaction channels.

#ifndef _PYNE_RXNAME_
#define _PYNE_RXNAME_
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <exception>
#include <stdlib.h>
#include <stdio.h>

#include "pyne.h"
#include "nucname.h"

/// Number of reactions supported by default.
#define NUM_RX_NAMES 540

namespace pyne
{
namespace rxname 
{
  extern std::string _names[NUM_RX_NAMES];  ///< Raw array of reaction names
  /// Set of reaction names, must be valid variable names.
  extern std::set<std::string> names;       
  /// Mapping from reaction ids to reaction names. 
  extern std::map<unsigned int, std::string> id_name;
  /// Mapping from reaction names to reaction ids. 
  extern std::map<std::string, unsigned int> name_id;
  /// Mapping between alternative names for reactions and the reaction id.
  extern std::map<std::string, unsigned int> altnames;
  /// Mapping from reaction ids to MT numbers.
  extern std::map<unsigned int, unsigned int> id_mt;
  /// Mapping from MT numbers to reaction names.
  extern std::map<unsigned int, unsigned int> mt_id;
  /// Mapping from reaction ids to labels (short descriptions).
  extern std::map<unsigned int, std::string> labels;
  /// Mapping from reaction ids to documentation strings (long descriptions).
  extern std::map<unsigned int, std::string> docs;
  /// Mapping from particle type, changes in Z number, and changes in A numbers 
  /// to reaction ids.  
  /// Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
  extern std::map<std::string, std::map<int, std::map<int, unsigned int> > > zadelta;
  /// A helper function to set the contents of the variables in this library.
  void * _fill_maps();  
  extern void * _;  ///< A dummy variable used when calling #_fill_maps().

  /// \name Hash Functions
  /// \{
  /// Custom hash function for reaction name to reaction ids.
  /// This functions will not return a value less than 1000, effectively reserving
  /// space for the MT numbers.
  unsigned int hash(std::string s);
  unsigned int hash(const char * s);
  /// \}

  /// \name Name Functions
  /// \{
  /// Returns the canonical name of a reaction channel.
  /// \param n Integer input of possible reaction, nominally an id or MT number.
  /// \param s String input of possible reaction, often a reaction or alternate name.
  /// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
  ///                 an integer it must be in id form.
  /// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
  ///               an integer it must be in id form.
  /// \param z Flag for incident particle type.
  ///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
  std::string name(int n);
  std::string name(unsigned int n);
  std::string name(char * s);
  std::string name(std::string s);
  std::string name(int from_nuc, int to_nuc, std::string z="n");
  std::string name(int from_nuc, std::string to_nuc, std::string z="n");
  std::string name(std::string from_nuc, int to_nuc, std::string z="n");
  std::string name(std::string from_nuc, std::string to_nuc, std::string z="n");
  /// \}

  /// \name ID Functions
  /// \{
  /// Returns the recation id of a reaction channel.  This id has been precomputed  
  /// from the hash of the name.
  /// \param x Input reaction specification, may be a reaction name, alternate name,
  ///          an id, or an MT number.
  /// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
  ///                 an integer it must be in id form.
  /// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
  ///               an integer it must be in id form.
  /// \param z Flag for incident particle type.
  ///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
  unsigned int id(int x);
  unsigned int id(unsigned int x);
  unsigned int id(char * x);
  unsigned int id(std::string x);
  unsigned int id(int from_nuc, int to_nuc, std::string z="n");
  unsigned int id(int from_nuc, std::string to_nuc, std::string z="n");
  unsigned int id(std::string from_nuc, int to_nuc, std::string z="n");
  unsigned int id(std::string from_nuc, std::string to_nuc, std::string z="n");
  /// \}

  /// \name MT Number Functions
  /// \{
  /// Returns the MT number of a reaction channel, if available.
  /// \param x Input reaction specification, may be a reaction name, alternate name,
  ///          an id, or an MT number.
  /// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
  ///                 an integer it must be in id form.
  /// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
  ///               an integer it must be in id form.
  /// \param z Flag for incident particle type.
  ///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
  unsigned int mt(int x);
  unsigned int mt(unsigned int x);
  unsigned int mt(char * x);
  unsigned int mt(std::string x);
  unsigned int mt(int from_nuc, int to_nuc, std::string z="n");
  unsigned int mt(int from_nuc, std::string to_nuc, std::string z="n");
  unsigned int mt(std::string from_nuc, int to_nuc, std::string z="n");
  unsigned int mt(std::string from_nuc, std::string to_nuc, std::string z="n");
  /// \}

  //// \name Label Functions
  /// \{
  /// Returns a short description of a reaction channel.
  /// \param x Input reaction specification, may be a reaction name, alternate name,
  ///          an id, or an MT number.
  /// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
  ///                 an integer it must be in id form.
  /// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
  ///               an integer it must be in id form.
  /// \param z Flag for incident particle type.
  ///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
  std::string label(int x);
  std::string label(unsigned int x);
  std::string label(char * x);
  std::string label(std::string x);
  std::string label(int from_nuc, int to_nuc, std::string z="n");
  std::string label(int from_nuc, std::string to_nuc, std::string z="n");
  std::string label(std::string from_nuc, int to_nuc, std::string z="n");
  std::string label(std::string from_nuc, std::string to_nuc, std::string z="n");
  /// \}

  /// \name Documentation Functions
  /// \{
  /// Returns a short description of a reaction channel.
  /// \param x Input reaction specification, may be a reaction name, alternate name,
  ///          an id, or an MT number.
  /// \param from_nuc Initial target nuclide prior to reaction.  When \a from_nuc is
  ///                 an integer it must be in id form.
  /// \param to_nuc Target nuclide after reaction occurs.  When \a to_nuc is
  ///               an integer it must be in id form.
  /// \param z Flag for incident particle type.
  ///          Particle flags are 'n', 'p', 'd', 't', 'He3', 'a', 'gamma', and 'decay'.
  std::string doc(int x);
  std::string doc(unsigned int x);
  std::string doc(char * x);
  std::string doc(std::string x);
  std::string doc(int from_nuc, int to_nuc, std::string z="n");
  std::string doc(int from_nuc, std::string to_nuc, std::string z="n");
  std::string doc(std::string from_nuc, int to_nuc, std::string z="n");
  std::string doc(std::string from_nuc, std::string to_nuc, std::string z="n");
  /// \}


  /// Custom exception for declaring a value not to be a valid reaction.  
  class NotAReaction : public std::exception
  {
  public:

    /// default constructor
    NotAReaction () {};

    /// default destructor
    ~NotAReaction () throw () {};

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    NotAReaction(std::string wasptr, std::string nowptr)
    {
       rxwas = wasptr;
       rxnow = nowptr;
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    NotAReaction(std::string wasptr, int nowptr)
    {
      rxwas = wasptr;
      rxnow = pyne::to_str(nowptr);
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    NotAReaction(int wasptr, std::string nowptr)
    {
      rxwas = pyne::to_str(wasptr);
      rxnow = nowptr;
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    NotAReaction(int wasptr, int nowptr)
    {
      rxwas = pyne::to_str(wasptr);
      rxnow = pyne::to_str(nowptr);
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    NotAReaction(std::string wasptr, unsigned int nowptr)
    {
      rxwas = wasptr;
      rxnow = pyne::to_str(nowptr);
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    NotAReaction(unsigned int wasptr, std::string nowptr)
    {
      rxwas = pyne::to_str(wasptr);
      rxnow = nowptr;
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    NotAReaction(unsigned int wasptr, unsigned int nowptr)
    {
      rxwas = pyne::to_str(wasptr);
      rxnow = pyne::to_str(nowptr);
    };

    /// Returns a helpful error message containing prior and current reaction state.
    virtual const char* what() const throw()
    {
      std::string narxstr ("Not a reaction! ");
      if (!rxwas.empty())
        narxstr += rxwas;

      if (!rxnow.empty())
      {
        narxstr += " --> "; 
        narxstr += rxnow;
      }
      return (const char *) narxstr.c_str();
    };

  private:
    std::string rxwas;  ///< previous reaction state
    std::string rxnow;  ///< current reaction state
  };



  /// Custom exception for declaring a value not to be of ambiquous reaction form.
  class IndeterminateReactionForm : public std::exception
  {
  public:

    /// default constructor
    IndeterminateReactionForm () {};

    /// default destructor
    ~IndeterminateReactionForm () throw () {};

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    IndeterminateReactionForm(std::string wasptr, std::string nowptr)
    {
       rxwas = wasptr;
       rxnow = nowptr;
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    IndeterminateReactionForm(std::string wasptr, int nowptr)
    {
      rxwas = wasptr;
      rxnow = pyne::to_str(nowptr);
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    IndeterminateReactionForm(int wasptr, std::string nowptr)
    {
      rxwas = pyne::to_str(wasptr);
      rxnow = nowptr;
    };

    /// Constructor using original reaction (\a wasptr) and the eventual state
    /// that PyNE calculated (\a nowptr).
    IndeterminateReactionForm(int wasptr, int nowptr)
    {
      rxwas = pyne::to_str(wasptr);
      rxnow = pyne::to_str(nowptr);
    };

    /// Returns a helpful error message containing prior and current reaction state.
    virtual const char* what() const throw()
    {
      std::string INFEstr ("Indeterminate reaction form: ");
      if (!rxwas.empty())
        INFEstr += rxwas;

      if (!rxnow.empty())
      {
        INFEstr += " --> "; 
        INFEstr += rxnow;
      }
      return (const char *) INFEstr.c_str();
    }

  private:
    std::string rxwas;  ///< previous reaction state
    std::string rxnow;  ///< current reaction state
  };
};
};

#endif
