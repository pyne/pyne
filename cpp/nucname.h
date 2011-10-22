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
  /*** Constructs the name to zz Dictionary ***/
  typedef std::string name_t;
  typedef int zz_t;

  typedef std::map<name_t, zz_t> name_zz_t; 
  typedef name_zz_t::iterator name_zz_iter; 
  name_zz_t get_name_zz();
  extern name_zz_t name_zz; 

  /*** Constructs zz to name dictionary **/
  typedef std::map<zz_t, name_t> zzname_t; 
  typedef zzname_t::iterator zzname_iter; 
  zzname_t get_zz_name();
  extern zzname_t zz_name;

  /******************************************/
  /*** Define useful elemental group sets ***/
  /******************************************/
  typedef std::set<name_t> name_group;
  typedef name_group::iterator name_group_iter;

  typedef std::set<zz_t> zz_group;
  typedef zz_group::iterator zz_group_iter;

  zz_group name_to_zz_group (name_group);

  extern name_t LAN_array[15];
  extern name_group LAN;
  extern zz_group lan;

  extern name_t ACT_array[15];
  extern name_group ACT;
  extern zz_group act;

  extern name_t TRU_array[19];
  extern name_group TRU;
  extern zz_group tru;

  extern name_t MA_array[10];
  extern name_group MA;
  extern zz_group ma;

  extern name_t FP_array[88];
  extern name_group FP;
  extern zz_group fp;


  /******************/
  /*** Exceptions ***/
  /******************/

  class NotANuclide : public std::exception
  {
  public:
    NotANuclide () {};
    ~NotANuclide () throw () {};

    NotANuclide(std::string wasptr, std::string nowptr)
    {
       nucwas = wasptr;
       nucnow = nowptr;
    };

    NotANuclide(std::string wasptr, int nowptr)
    {
      nucwas = wasptr;
      nucnow = pyne::to_str(nowptr);
    };

    NotANuclide(int wasptr, std::string nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = nowptr;
    };

    NotANuclide(int wasptr, int nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = pyne::to_str(nowptr);
    };

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
    std::string nucwas;
    std::string nucnow;
  };



  class IndeterminateNuclideForm : public std::exception
  {
  public:
    IndeterminateNuclideForm () {};
    ~IndeterminateNuclideForm () throw () {};

    IndeterminateNuclideForm(std::string wasptr, std::string nowptr)
    {
       nucwas = wasptr;
       nucnow = nowptr;
    };

    IndeterminateNuclideForm(std::string wasptr, int nowptr)
    {
      nucwas = wasptr;
      nucnow = pyne::to_str(nowptr);
    };

    IndeterminateNuclideForm(int wasptr, std::string nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = nowptr;
    };

    IndeterminateNuclideForm(int wasptr, int nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = pyne::to_str(nowptr);
    };

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
    std::string nucwas;
    std::string nucnow;
  };


  /********************/
  /*** Current Form ***/
  /********************/
  std::string current_form(std::string);
  std::string current_form(int);


  /************************/
  /*** zzaaam functions ***/
  /************************/
  int zzaaam(int);
  int zzaaam(char *);
  int zzaaam(std::string);


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
