// Material.h
// Class that defines material flows for the Fuel Cycle.  Input into most FCComp objects.
// Effectively an normalized isotopic linked list with associated mass.
// Contains other functions for mixing mass streams and setting up sub streams.
//
// Initial Author --- Anthony Scopatz

#if !defined(_PYNE_MATERIAL_)
#define _PYNE_MATERIAL_

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>

#define JSON_IS_AMALGAMATION
#include <json/json-forwards.h>
#include <json/json.h>
#include "h5wrap.h"

#include "pyne.h"
#include "nucname.h"
#include "data.h"

namespace pyne
{
  // Set Type Definitions
  typedef std::map<int, double> comp_map;
  typedef comp_map::iterator comp_iter;


  class Material
  {
  // Parent Class for Handling Mass Streams"
  protected:
    //Material function definitions
    double get_comp_sum ();

  public:
    void norm_comp ();

    // Persistence functions.
    void _load_comp_protocol0(hid_t, std::string, int);
    void _load_comp_protocol1(hid_t, std::string, int);

    void from_hdf5(char *, char *, int=-1, int=1);
    void from_hdf5(std::string, std::string="/material", int=-1, int=1);

    void write_hdf5(char *, char *, char *, float=-0.0, int=100);
    void write_hdf5(std::string, std::string="/material", std::string nuclist="/nuc_zz", float=-0.0, int=100);

    void from_text(char *);
    void from_text(std::string);

    void write_text(char *);
    void write_text(std::string);

    //Fundemental mass stream data
    comp_map comp;
    double mass;
    double density;
    double atoms_per_mol;
    Json::Value attrs;

    //Material Constructors
    Material ();
    Material(comp_map, double=-1.0, double=-1.0, double=-1.0,
             Json::Value=Json::Value(Json::objectValue));
    Material(char *, double=-1.0, double=-1.0, double=-1.0,
             Json::Value=Json::Value(Json::objectValue));
    Material(std::string, double=-1.0, double=-1.0, double=-1.0,
             Json::Value=Json::Value(Json::objectValue));
    ~Material ();

    //Material function definitions
    void normalize ();
    comp_map mult_by_mass();
    double molecular_weight(double=-1.0);

    //Sub-Stream Computation
    Material sub_mat(std::set<int>);
    Material sub_mat(std::set<std::string>);

    Material set_mat(std::set<int>, double);
    Material set_mat(std::set<std::string>, double);

    Material del_mat(std::set<int>);
    Material del_mat(std::set<std::string>);

    Material sub_range(int=0, int=10000000);
    Material set_range(int=0, int=10000000, double=0.0);
    Material del_range(int=0, int=10000000);

    Material sub_u();
    Material sub_pu();
    Material sub_lan();
    Material sub_act();
    Material sub_tru();
    Material sub_ma();
    Material sub_fp();

    // Atom fraction functions
    std::map<int, double> to_atom_frac();
    void from_atom_frac(std::map<int, double>);

    //Overloaded Operators
    Material operator+ (double);
    Material operator+ (Material);
    Material operator* (double);
    Material operator/ (double);
  };

  std::ostream& operator<< (std::ostream& os, Material mat);

  typedef struct material_struct {
    double mass;
    double density;
    double atoms_per_mol;
    double comp [];
  } material_struct;


  /******************/
  /*** Exceptions ***/
  /******************/

  class MaterialProtocolError: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Invalid loading protocol number; please use 0 or 1.";
    };
  };

// End pyne namespace
};

#endif
