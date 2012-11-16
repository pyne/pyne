// Converts between naming conventions for reaction channels.

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

namespace pyne
{
namespace rxname 
{
  extern std::string _names[15];
  extern std::set<std::string> names;

  /******************/
  /*** Exceptions ***/
  /******************/
  class NotAReaction : public std::exception
  {
  public:
    NotAReaction () {};
    ~NotAReaction () throw () {};

    NotAReaction(std::string wasptr, std::string nowptr)
    {
       nucwas = wasptr;
       nucnow = nowptr;
    };

    NotAReaction(std::string wasptr, int nowptr)
    {
      nucwas = wasptr;
      nucnow = pyne::to_str(nowptr);
    };

    NotAReaction(int wasptr, std::string nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = nowptr;
    };

    NotAReaction(int wasptr, int nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = pyne::to_str(nowptr);
    };

    virtual const char* what() const throw()
    {
      std::string narxstr ("Not a reaction! ");
      if (!nucwas.empty())
        narxstr += nucwas;

      if (!nucnow.empty())
      {
        narxstr += " --> "; 
        narxstr += nucnow;
      }
      return (const char *) narxstr.c_str();
    };

  private:
    std::string nucwas;
    std::string nucnow;
  };



  class IndeterminateReactionForm : public std::exception
  {
  public:
    IndeterminateReactionForm () {};
    ~IndeterminateReactionForm () throw () {};

    IndeterminateReactionForm(std::string wasptr, std::string nowptr)
    {
       nucwas = wasptr;
       nucnow = nowptr;
    };

    IndeterminateReactionForm(std::string wasptr, int nowptr)
    {
      nucwas = wasptr;
      nucnow = pyne::to_str(nowptr);
    };

    IndeterminateReactionForm(int wasptr, std::string nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = nowptr;
    };

    IndeterminateReactionForm(int wasptr, int nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = pyne::to_str(nowptr);
    };

    virtual const char* what() const throw()
    {
      std::string INFEstr ("Indeterminate reaction form: ");
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


  /************************/
  /*** name functions *****/
  /************************/
  std::string name(int);
  std::string name(char *);
  std::string name(std::string);

  /********************/
  /*** mt functions ***/
  /********************/
/*
  int mt(int);
  int mt(char *);
  int mt(std::string);
*/

};
};

#endif
