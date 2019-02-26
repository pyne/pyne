#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>

namespace pyne {
class MaterialLibrary() {
 protected:
  // The actual library
  std::map<std::string, pyne::Material> material_library;  // material library

 private:
  // turns the filename string into the full file path
  std::string get_full_filepath(char* filename);
  // turns the filename string into the full file path
  std::string get_full_filepath(std::string filename);
  // make sure the file, filename exists
  bool check_file_exists(std::string filename);

 public:
  // materialLibrary constructor
  MaterialLibrary();  //< empty constryctor

  /// Constructor from file
  /// \param filename path to file on disk, this file may be either in plaintext
  ///                 or HDF5 format.
  MaterialLibrary(string filename, std::string datapath = "/materials");

  ~MaterialLibrary();  //< default destructor

  /**
   * \brief loads the pyne materials in map of name vs Material
   * \param[in] filename of the h5m file
   * \return std::map of material name vs Material object
   */
  void load_materials(std::string filename,
                      std::string datapath = "/materials");

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

}  // end MaterialLibrary class header
}  // end of pyne namespace
