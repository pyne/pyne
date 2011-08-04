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

#include "H5Cpp.h"
#include "h5wrap.h"

#include "pyne.h"
#include "nucname.h"

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

    void load_from_hdf5 (char *, char *, int = -1);
    void load_from_hdf5 (std::string, std::string="/material", int = -1);

    void load_from_text (char *);
    void load_from_text (std::string);

    //Fundemental mass stream data
    comp_map comp;
    double mass;
    std::string name;

    //Material Constructors
    Material ();
    Material (comp_map, double = -1.0, std::string = "");
    Material (char *, double = -1.0, std::string = "");
    Material (std::string, double = -1.0, std::string = "");
    ~Material ();

    //Material function definitions
    void normalize ();
    comp_map mult_by_mass();
    double atomic_weight();

    //Sub-Stream Computation
    Material sub_mat (std::set<int>,  std::string n = "");
    Material sub_mat (std::set<std::string>,  std::string n = "");
    Material sub_u   (std::string = "");
    Material sub_pu  (std::string = "");
    Material sub_lan (std::string = "");
    Material sub_act (std::string = "");
    Material sub_tru (std::string = "");
    Material sub_ma  (std::string = "");
    Material sub_fp  (std::string = "");

    //Overloaded Operators
    Material operator+ (double);
    Material operator+ (Material);
    Material operator* (double);
    Material operator/ (double);
  };

std::ostream& operator<< (std::ostream& os, Material mat);

// End pyne namespace
};

#endif
