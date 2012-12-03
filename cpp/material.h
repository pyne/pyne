/// \file material.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief The ever-important material class and related helpers.
/// 
/// The material class is effectively a normalized nuclide linked list with 
/// associated mass, density, atoms per mol, and metadata.  However, this 
/// implementation also contains other functions for mixing materials and generating 
/// related materials.

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
  typedef std::map<int, double> comp_map; ///< Nuclide-mass composition map type
  typedef comp_map::iterator comp_iter;   ///< Nuclide-mass composition iter type

  /// Material composed of nuclides.
  class Material
  {
  protected:

    /// Computes the total mass stored in the composition.
    double get_comp_sum ();

  public:

    /// Normalizes the mass values in the composition.
    void norm_comp ();

    // Persistence functions.

    /// Loads the matrial composition from an HDF5 file according to the layout 
    /// defined by protocol 0.  This protocol is depratacted.
    /// \param db HDF5 id for the open HDF5 file.
    /// \param datapath Path to the base node for the material in \a db.
    /// \param row The index to read out, may be negative.
    void _load_comp_protocol0(hid_t db, std::string datapath, int row);

    /// Loads the matrial composition from an HDF5 file according to the layout 
    /// defined by protocol 1.  This protocol should be used in favor of protocol 0.
    /// \param db HDF5 id for the open HDF5 file.
    /// \param datapath Path to the base node for the material in \a db.
    /// \param row The index to read out, may be negative.
    void _load_comp_protocol1(hid_t db, std::string datapath, int row);

    /// Loads a material from an HDF5 file into this object.
    /// \param filename Path on disk to the HDF5 file.
    /// \param datapath Path to the the material in the file.
    /// \param row The index to read out, may be negative.
    /// \param protocol Flag for layout of material on disk.
    void from_hdf5(char * filename, char * datapath, int row=-1, int protocol=1);

    /// Loads a material from an HDF5 file into this object.
    /// \param filename Path on disk to the HDF5 file.
    /// \param datapath Path to the the material in the file.
    /// \param row The index to read out, may be negative.
    /// \param protocol Flag for layout of material on disk.
    void from_hdf5(std::string filename, std::string datapath="/material", 
                                                          int row=-1, int protocol=1);

    /// Writes this material out to an HDF5 file.  
    /// This happens according to protocol 1.
    /// \param filename Path on disk to the HDF5 file.
    /// \param datapath Path to the the material in the file.
    /// \param nucpath Path to the nuclides set in the file.
    /// \param row The index to read out, may be negative. Also note that this is a
    ///            float.  A value of -0.0 indicates that the material should be 
    ///            appended to the end of the dataset.
    /// \param chunksize The chunksize for all material data on disk.
    void write_hdf5(char * filename, char * datapath, char * nucpath, float row=-0.0, 
                                                                    int chunksize=100);
    /// Writes this material out to an HDF5 file.  
    /// This happens according to protocol 1.
    /// \param filename Path on disk to the HDF5 file.
    /// \param datapath Path to the the material in the file.
    /// \param nucpath Path to the nuclides set in the file.
    /// \param row The index to read out, may be negative. Also note that this is a
    ///            float.  A value of -0.0 indicates that the material should be 
    ///            appended to the end of the dataset.
    /// \param chunksize The chunksize for all material data on disk.
    void write_hdf5(std::string filename, std::string datapath="/material", 
                    std::string nucpath="/nuc_zz", float row=-0.0, int chunksize=100);

    /// Reads data from a plaintext file at \a filename into this Material instance.
    void from_text(char * filename);
    /// Reads data from a plaintext file at \a filename into this Material instance.
    void from_text(std::string filname);

    /// Writes the Material out to a simple plaintext file readable by from_text().
    void write_text(char * filename);
    /// Writes the Material out to a simple plaintext file readable by from_text().
    void write_text(std::string filename);

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
    Material expand_elements();

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
