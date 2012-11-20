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
#include "nucname.h"

#define NUM_RX_NAMES 328

namespace pyne
{
namespace rxname 
{
  extern std::string _names[NUM_RX_NAMES];
  extern std::set<std::string> names;
  extern std::map<std::string, unsigned int> altnames;
  extern std::map<unsigned int, std::string> id_name;
  extern std::map<std::string, unsigned int> name_id;
  extern std::map<unsigned int, unsigned int> id_mt;
  extern std::map<unsigned int, unsigned int> mt_id;
  extern std::map<unsigned int, std::string> labels;
  extern std::map<unsigned int, std::string> docs;
  extern std::map<std::string, std::map<int, std::map<int, unsigned int> > > zadelta;
  void * _fill_maps();
  extern void * _;

  unsigned int hash(std::string);
  unsigned int hash(const char *);

  /************************/
  /*** name functions *****/
  /************************/
  std::string name(int);
  std::string name(unsigned int);
  std::string name(char *);
  std::string name(std::string);
  std::string name(int, int, std::string="n");
  std::string name(int, std::string, std::string="n");
  std::string name(std::string, int, std::string="n");
  std::string name(std::string, std::string, std::string="n");

  /**********************/
  /*** id functions *****/
  /**********************/
  unsigned int id(int);
  unsigned int id(unsigned int);
  unsigned int id(char *);
  unsigned int id(std::string);
  unsigned int id(int, int, std::string="n");
  unsigned int id(int, std::string, std::string="n");
  unsigned int id(std::string, int, std::string="n");
  unsigned int id(std::string, std::string, std::string="n");

  /********************/
  /*** MT functions ***/
  /********************/
  unsigned int mt(int);
  unsigned int mt(unsigned int);
  unsigned int mt(char *);
  unsigned int mt(std::string);
  unsigned int mt(int, int, std::string="n");
  unsigned int mt(int, std::string, std::string="n");
  unsigned int mt(std::string, int, std::string="n");
  unsigned int mt(std::string, std::string, std::string="n");

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

    NotAReaction(std::string wasptr, unsigned int nowptr)
    {
      nucwas = wasptr;
      nucnow = pyne::to_str(nowptr);
    };

    NotAReaction(unsigned int wasptr, std::string nowptr)
    {
      nucwas = pyne::to_str(wasptr);
      nucnow = nowptr;
    };

    NotAReaction(unsigned int wasptr, unsigned int nowptr)
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
};
};

#endif
