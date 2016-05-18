#include "../pyne/pyne.h"
#include <stdlib.h>
#include <map>
#include <string>

//===========================================================================//
/**
 * \class UWUW
 * \brief Defines the UWUW interface
 *
 * UWUW is a base class that defines the variables and methods
 * needed to load University of Wisconsin Unified Workflow data from the filename
 * pointed to. It contains the following private functions,
 *
 *    std::string get_full_filepath(char *filename);
 *    std::string get_full_filepath(std::string filename);
 *    std::map<std::string, pyne::Material> load_pyne_materials(std::string filename);
 *    std::map<std::string, pyne::Material> load_pyne_tallies(std::string filename);
 *
 * Note there are two constructors, one that takes std::string filename as argument and
 * the other char *, so both forms of strings are handled.
 *
 * The filename must be marked up appropriately with the material objects that are to be
 * used in the workflow, currently that is performed using a Python script, uwuw_preproc
 * to insert these. We use the python write_hdf5 option, with the hdf5 filepath to be
 * "/materials", so when load_hdf5() is called in this class, it looks for material objects
 * in the "/materials" path.
 *
 * Tallies are saved to "/tally", so the load_hdf5() method also looks in that path "/tally"
 * otherwise it will return a map of 0 size.
 */
//===========================================================================//

class UWUW
{
 public:
  /**
   * \brief Constructor
   * Empty constructor
   */
  UWUW(); // empty constructor

  /**
   * \brief Instantiates instance of the UWUW class, populates 3 member variables,
   * full_filepath, tally_library and material_library.
   *
   * \param[in] filename, the filename containing pyne material & tally objects
   *
   * the get_full_filepath() method populates a string containing the filename passed in
   * with the the path of the file appended to it.
   * the load_pyne_materials() methods populates a map of material name vs Material object
   * the load_pyne_tallies() methods populates a map of tally name vs tally object
   */
  UWUW(char* filename); // c style filename

  /**
   * \brief Instantiates instance of the UWUW class, populates 3 member variables,
   * full_filepath, tally_library and material_library.
   *
   * \param[in] filename, the filename containing pyne material & tally objects
   *
   * the get_full_filepath() method populates a string containing the filename passed in
   * with the the path of the file appended to it.
   * the load_pyne_materials() methods populates a map of material name vs Material object
   * the load_pyne_tallies() methods populates a map of tally name vs tally object
   */
  UWUW(std::string filename); // normal constructor

  /**
   * \brief Standard destructor
   */
  ~UWUW(); // destructor

  // Member variables

  /**
   * \brief material_library is a std::map of <material name, Material object>
   * iterate through this map to get its contents
   */
  std::map<std::string, pyne::Material> material_library; // material library
  /**
   * \brief tally_library is a std::map of <tally name, Tally object>
   * iterate through this map to get its contents
   */
  std::map<std::string, pyne::Tally> tally_library; // tally library
  /**
   * \brief full filepath is variable set to /path/+filename
   * This is needed to set the pyne::nucdata path
   */
  std::string full_filepath;
  /**
   * \brief the number of members in the tally library
   */
  int num_tallies;
  /**
   * \brief the number of members in the material libary
   */
  int num_materials;


 private:
  // turns the filename string into the full file path
  std::string get_full_filepath(char *filename);
  // turns the filename string into the full file path
  std::string get_full_filepath(std::string filename);
  // make sure the file, filename exists
  bool check_file_exists(std::string filename);

  /**
   * \brief loads the pyne materials in map of name vs Material
   * \param[in] filename of the h5m file
   * \return std::map of material name vs Material object
   */
 public:
  std::map<std::string, pyne::Material> load_pyne_materials(std::string filename, std::string datapath = "/materials");

  /**
   * \brief loads the pyne tallies in map of name vs Material
   * \param[in] filename of the h5m file
   * \return std::map of tally name vs Tally object
   */
  std::map<std::string, pyne::Tally> load_pyne_tallies(std::string filename, std::string datapath = "/tally");

  /**
   * \brief determines in that datapath exists in the hdf5 file
   * \param[in] filename of the h5m file
   * \param[in] the datapath of we would like to test
   * \return true/false
   */
  bool hdf5_path_exists(std::string filename, std::string datapath);

  /**
   * \brief determines the length of an hdf5 data table
   * \param[in] filename of the h5m file
   * \param[in] the datapath of we would like to read
   * \return the number of elements to the array
   */
  int get_length_of_table(std::string filename, std::string datapath);

};
