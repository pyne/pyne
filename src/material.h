/// \brief The ever-important material class and related helpers.
///
/// The material class is effectively a normalized nuclide linked list with
/// associated mass, density, atoms per mol, and metadata.  However, this
/// implementation also contains other functions for mixing materials and generating
/// related materials.

#ifndef PYNE_MR34UE5INRGMZK2QYRDWICFHVM
#define PYNE_MR34UE5INRGMZK2QYRDWICFHVM

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>	// std::ostringstream

#if !defined(JSON_IS_AMALGAMATION)
  #define JSON_IS_AMALGAMATION
#endif

#ifndef PYNE_IS_AMALGAMATED
#include "json-forwards.h"
#include "json.h"
#include "h5wrap.h"
#include "utils.h"
#include "nucname.h"
#include "data.h"
#endif

namespace pyne
{
  // Set Type Definitions
  typedef std::map<int, double> comp_map; ///< Nuclide-mass composition map type
  typedef comp_map::iterator comp_iter;   ///< Nuclide-mass composition iter type

  // These 37 strings are predefined FLUKA materials. 
  // Materials not on this list requires a MATERIAL card. 
  static std::string fluka_mat_strings[] = {
   "BLCKHOLE", "VACUUM",   "HYDROGEN", "HELIUM",   "BERYLLIU", "CARBON", 
   "NITROGEN", "OXYGEN",   "MAGNESIU", "ALUMINUM", "IRON",     "COPPER", 
   "SILVER",   "SILICON",  "GOLD",     "MERCURY",  "LEAD",     "TANTALUM", 
   "SODIUM",   "ARGON",    "CALCIUM",  "TIN",      "TUNGSTEN", "TITANIUM", 
   "NICKEL",   "WATER",    "POLYSTYR", "PLASCINT", "PMMA",     "BONECOMP", 
   "BONECORT", "MUSCLESK", "MUSCLEST", "ADTISSUE", "KAPTON", "POLYETHY", "AIR"
  };
  static int FLUKA_MAT_NUM = 37;

  /// Material composed of nuclides.
  class Material
  {
  protected:

    /// Computes the total mass stored in the composition.
    double get_comp_sum ();

  public:

    // Material Constructors
    Material ();  ///< empty constructor
    /// Constructor from composition map
    /// \param cm composition map
    /// \param m mass value, the mass is set to the sum of the values in the
    ///          composition if \a m is negative.
    /// \param d density value
    /// \param apm atoms per mole
    /// \param attributes initial metadata
    Material(comp_map cm, double m=-1.0, double d=-1.0, double apm=-1.0,
             Json::Value attributes=Json::Value(Json::objectValue));
    /// Constructor from file
    /// \param filename path to file on disk, this file may be either in plaintext
    ///                 or HDF5 format.
    /// \param m mass value, the mass is set to the sum of the values in the
    ///          composition if \a m is negative,
    ///          may be overridden by the value from disk.
    /// \param d density value,
    ///          may be overridden by the value from disk.
    /// \param apm atoms per mole,
    ///          may be overridden by the value from disk.
    /// \param attributes initial metadata,
    ///          may be overridden by the value from disk.
    Material(char * filename, double m=-1.0, double d=-1.0, double apm=-1.0,
             Json::Value attributes=Json::Value(Json::objectValue));
    /// Constructor from file
    /// \param filename path to file on disk, this file may be either in plaintext
    ///                 or HDF5 format.
    /// \param m mass value, the mass is set to the sum of the values in the
    ///          composition if \a m is negative,
    ///          may be overridden by the value from disk.
    /// \param d density value,
    ///          may be overridden by the value from disk.
    /// \param apm atoms per mole,
    ///          may be overridden by the value from disk.
    /// \param attributes initial metadata,
    ///          may be overridden by the value from disk.
    Material(std::string filename, double m=-1.0, double d=-1.0, double apm=-1.0,
             Json::Value attributes=Json::Value(Json::objectValue));
    ~Material (); ///< default destructor

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
                    std::string nucpath="/nucid", float row=-0.0, int chunksize=100);

    /// Return an mcnp input deck record as a string
    std::string mcnp(std::string frac_type = "mass");
    /// 
    /// Return a fluka input deck MATERIAL card as a string
    std::string fluka(int id, std::string frac_type = "mass");
    /// Convenience function to tell whether a given name needs a material card
    bool not_fluka_builtin(std::string fluka_name);
    /// High level call to get details and call material_component(..)
    std::string fluka_material_str(int id);
    /// Intermediate level call to prepare final info and call material_line(..)
    std::string fluka_material_component(int fid, int nucid, 
                                         std::string fluka_name);
    /// Format information into a FLUKA material card
    std::string fluka_material_line(int znum, double atomic_mass, 
                              int fid, std::string fluka_name);
    /// Convenience function to format a single fluka field
    std::string fluka_format_field(float field);
    /// Return FLUKA compound card and the material card for the named compound
    /// but not the material cards of the components
    std::string fluka_compound_str(int id, std::string frac_type = "mass");

    /// Reads data from a plaintext file at \a filename into this Material instance.
    void from_text(char * filename);
    /// Reads data from a plaintext file at \a filename into this Material instance.
    void from_text(std::string filename);

    /// Writes the Material out to a simple plaintext file readable by from_text().
    void write_text(char * filename);
    /// Writes the Material out to a simple plaintext file readable by from_text().
    void write_text(std::string filename);

    /// Loads a JSON instance tree into this Material.
    void load_json(Json::Value);
    /// Dumps the Material out to a JSON instance tree.
    Json::Value dump_json();
    /// Reads data from a JSON file at \a filename into this Material instance.
    void from_json(char * filename);
    /// Reads data from a JSON file at \a filename into this Material instance.
    void from_json(std::string filname);
    /// Writes the Material out to a JSON file
    void write_json(char * filename);
    /// Writes the Material out to a JSON file
    void write_json(std::string filename);

    // Fundemental mass stream data
    /// composition, maps nuclides in id form to normalized mass weights.
    comp_map comp;
    double mass;  ///< mass (in arbitrary units) of the Material.
    double density; ///< density (in arbitrary units) of the Material.
    double atoms_per_molecule; ///< The number of atoms per molecule.
    /// container for arbitrary metadata, following the JSON rules.
    Json::Value metadata;

    // Material function definitions
    void normalize ();  ///< Normalizes the mass.
    /// Returns a composition map that has been unnormalized by multiplying each
    /// mass weight by the actual mass of the material.
    comp_map mult_by_mass();
    /// Calculates the atomic weight of this material based on the composition
    /// and the number of atoms per mol.  If \a apm is non-negative then it is
    /// used (and stored on the instance) as the atoms_per_molecule for this calculation.
    /// If \a apm and atoms_per_molecule on this instance are both negative, then the best
    /// guess value calculated from the normailized composition is used here.
    double molecular_mass(double apm=-1.0);
    /// Returns a copy of the current material where all natural elements in the
    /// composition are expanded to their natural isotopic abundances.
    Material expand_elements();
    // Returns a copy of the current material where all the isotopes of the elements
    // are added up, atomic-fraction-wise, unless they are in the exception set
    Material collapse_elements(std::set<int> exception_znum);
    // Wrapped version to facilitate calling from python
    Material collapse_elements(int **int_ptr_arry);
    // void print_material( pyne::Material test_mat);
    /// Computes, sets, and returns the mass density when \a num_dens is greater
    /// than or equal zero.  If \a num_dens is negative, this simply returns the
    /// current value of the density member variable.  You may also use / set the
    /// atoms per molecule (atoms_per_molecule) in this function using \a apm.
    double mass_density(double num_dens=-1.0, double apm=-1.0);
    /// Computes and returns the number density of the material using the
    /// mass density if \a mass_dens is greater than or equal to zero.  If
    /// \a mass_dens is negative, the denisty member variable is used instead.
    /// You may also use / set the atoms per molecule (atoms_per_molecule) in this
    /// function using \a apm.
    double number_density(double mass_dens=-1.0, double apm=-1.0);

    // Sub-Stream Computation
    /// Creates a sub-Material with only the nuclides present in \a nucset.
    /// Elements of this set may be either in id form or simple Z numbers.
    Material sub_mat(std::set<int> nucset);
    /// Creates a sub-Material with only the nuclides present in \a nucset.
    /// Elements of this set may be in any form.
    Material sub_mat(std::set<std::string> nucset);

    /// Creates a new Material with the mass weights for all nuclides in \a nucset
    /// set to \a value.
    /// Elements of \a nucset may be either in id form or simple Z numbers.
    Material set_mat(std::set<int> nucset, double value);
    /// Creates a new Material with the mass weights for all nuclides in \a nucset
    /// set to \a value.  Elements of \a nucset may be in any form.
    Material set_mat(std::set<std::string> nucset, double value);

    /// Creates a new Material with the all nuclides in \a nucset removed.
    /// Elements of \a nucset may be either in id form or simple Z numbers.
    Material del_mat(std::set<int> nucset);
    /// Creates a new Material with the all nuclides in \a nucset removed.
    /// Elements of \a nucset may be in any form.
    Material del_mat(std::set<std::string> nucset);

    /// Creates a sub-Material based on a range of id-form integers.
    Material sub_range(int lower=0, int upper=10000000);
    /// Creates a new Material with the mass weights for all nuclides in the id
    /// range set to \a value.
    Material set_range(int lower=0, int upper=10000000, double value=0.0);
    /// Creates a new Material with the all nuclides in the id range removed.
    Material del_range(int lower=0, int upper=10000000);

    /// Creates a sub-Material of only the given element. Assumes element is
    /// id form.
    Material sub_elem(int element);
    /// Creates a sub-Material of only lanthanides.
    Material sub_lan();
    /// Creates a sub-Material of only actinides.
    Material sub_act();
    /// Creates a sub-Material of only transuranics.
    Material sub_tru();
    /// Creates a sub-Material of only minor actinides.
    Material sub_ma();
    /// Creates a sub-Material of only fission products.
    Material sub_fp();

    // Atom fraction functions
    /// Returns a mapping of the nuclides in this material to their atom fractions.
    /// This calculation is based off of the materials molecular weight.
    std::map<int, double> to_atom_frac();
    /// Sets the composition, mass, and atoms_per_molecule of this material to those
    /// calculated from \a atom_fracs, a mapping of nuclides to atom fractions values.
    void from_atom_frac(std::map<int, double> atom_fracs);

    // Radioactive Material functions
    /// Returns a list of gamma-rays energies in keV and intensities in
    /// decays/s/atom material unnormalized
    std::vector<std::pair<double, double> > gammas();
    /// Returns a list of x-rays average energies in keV and intensities in
    /// decays/s material unnormalized
    std::vector<std::pair<double, double> > xrays();
    /// Returns a list of photon energies in keV and intensities in
    /// decays/s/atom material unnormalized
    std::vector<std::pair<double, double> > photons(bool norm);
    /// Takes a list of photon energies and intensities and normalizes them
    /// so the sum of the intensities is one
    std::vector<std::pair<double, double> > normalize_radioactivity(
      std::vector<std::pair<double, double> > unnormed);


    // Overloaded Operators
    /// Adds mass to a material instance.
    Material operator+ (double);
    /// Adds two materials together.
    Material operator+ (Material);
    /// Multiplies a material's mass.
    Material operator* (double);
    /// Divides a material's mass.
    Material operator/ (double);
  };

  /// Converts a Material to a string stream representation for canonical writing.
  /// This operator is also defined on inheritors of std::ostream
  std::ostream& operator<< (std::ostream& os, Material mat);

  /// A stuct for reprensenting fundemental data in a material.
  /// Useful for HDF5 representations.
  typedef struct material_data {
    double mass;  ///< material mass
    double density; ///< material density
    double atoms_per_mol; ///< material atoms per mole
    double comp []; ///< array of material composition mass weights.
  } material_data;

  /// Custom exception for invalid HDF5 protocol numbers
  class MaterialProtocolError: public std::exception
  {
    /// marginally helpful error message.
    virtual const char* what() const throw()
    {
      return "Invalid loading protocol number; please use 0 or 1.";
    };
  };

// End pyne namespace
};

#endif  // PYNE_MR34UE5INRGMZK2QYRDWICFHVM
