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

    void from_hdf5 (char *, char *, int = -1);
    void from_hdf5 (std::string, std::string="/material", int = -1);

    void write_hdf5 (char *, char *, char *, int = -1);
    void write_hdf5 (std::string, std::string="/material", std::string nuclist="/nuc_zz", int = -1);

    void from_text (char *);
    void from_text (std::string);

    //Fundemental mass stream data
    comp_map comp;
    double mass;
    std::string name;
    double atoms_per_mol;

    //Material Constructors
    Material ();
    Material (comp_map, double=-1.0, std::string="", double=-1.0);
    Material (char *, double=-1.0, std::string="", double=-1.0);
    Material (std::string, double=-1.0, std::string="", double=-1.0);
    ~Material ();

    //Material function definitions
    void normalize ();
    comp_map mult_by_mass();
    double molecular_weight(double=-1.0);

    //Sub-Stream Computation
    Material sub_mat(std::set<int>, std::string n = "");
    Material sub_mat(std::set<std::string>, std::string n = "");

    Material set_mat(std::set<int>, double, std::string n = "");
    Material set_mat(std::set<std::string>, double, std::string n = "");

    Material del_mat(std::set<int>,  std::string n = "");
    Material del_mat(std::set<std::string>,  std::string n = "");

    Material sub_range(int=0, int=10000000, std::string n = "");
    Material set_range(int=0, int=10000000, double=0.0, std::string n = "");
    Material del_range(int=0, int=10000000, std::string n = "");

    Material sub_u   (std::string = "");
    Material sub_pu  (std::string = "");
    Material sub_lan (std::string = "");
    Material sub_act (std::string = "");
    Material sub_tru (std::string = "");
    Material sub_ma  (std::string = "");
    Material sub_fp  (std::string = "");

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
    char name [20];
    double mass;
    double atoms_per_mol;
    double comp [];
  } material_struct;

// End pyne namespace
};

#endif
