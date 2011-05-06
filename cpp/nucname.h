// Converts between naming conventions for nuclides.
// zzaaam is for numerals only (923350).
// LLAAAM is for letters  as well (U-235).
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

namespace nucname 
{
  /*** Constructs the LL to zz Dictionary ***/
  typedef std::string LL_t;
  typedef int zz_t;

  typedef std::map<LL_t, zz_t> LLzz_t; 
  typedef LLzz_t::iterator LLzz_iter; 
  LLzz_t get_LLzz();
  extern LLzz_t LLzz; 

  /*** Constructs zz to LL dictionary **/
  typedef std::map<zz_t, LL_t> zzLL_t; 
  typedef zzLL_t::iterator zzLL_iter; 
  zzLL_t get_zzLL();
  extern zzLL_t zzLL;

  /******************************************/
  /*** Define useful elemental group sets ***/
  /******************************************/
  typedef std::set<LL_t> LL_group;
  typedef LL_group::iterator LL_group_iter;

  typedef std::set<zz_t> zz_group;
  typedef zz_group::iterator zz_group_iter;

  zz_group LL_to_zz_group (LL_group);

  extern LL_t LAN_array[15];
  extern LL_group LAN;
  extern zz_group lan;

  extern LL_t ACT_array[23];
  extern LL_group ACT;
  extern zz_group act;

  extern LL_t TRU_array[19];
  extern LL_group TRU;
  extern zz_group tru;

  extern LL_t MA_array[18];
  extern LL_group MA;
  extern zz_group ma;

  LL_group get_FP();
  extern LL_group FP;
  LL_t * get_FP_array(LL_group);
  extern LL_t * FP_array;
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
  /*** LLAAAM functions ***/
  /************************/
  std::string LLAAAM(int);
  std::string LLAAAM(char *);
  std::string LLAAAM(std::string);


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




  /************************/
  /*** Helper Functions ***/
  /************************/
  double nuc_weight(int);
  double nuc_weight(char *);
  double nuc_weight(std::string);
}

#endif
