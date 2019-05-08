#ifndef PYNE_MATERIAL_LIBRARY
#define PYNE_MATERIAL_LIBRARY

#include <stdlib.h>
#include <iostream>
#include <unordered_map>
#include <string>
#include <string>
#include <vector>

#if !defined(JSON_IS_AMALGAMATION)
#define JSON_IS_AMALGAMATION
#endif

#ifndef PYNE_IS_AMALGAMATED
#include "material.h"
#endif

namespace pyne {

typedef std::unordered_map<std::string, pyne::Material> mat_map;
typedef std::set<std::string> matname_set;
typedef std::set<int> nuc_set;

class MaterialLibrary {
 protected:
  // The actual library
  mat_map material_library;  // material library
  matname_set matlist;
  nuc_set nuclist;


 private:
  /**
   * \brief Add the nuclide form the material to the list of nuclides
   * \param material from which add the nuclide to the list
  */
  void append_to_nuclist(pyne::Material mat);
  /**
   * \brief determines in that datapath exists in the hdf5 file
   * \param[in] filename of the h5m file
   * \param[in] the datapath of we would like to test
   * \return true/false
   */
  bool hdf5_path_exists(const std::string& filename,
                        const std::string& datapath);

  /**
   * \brief determines the length of an hdf5 data table
   * \param[in] filename of the h5m file
   * \param[in] the datapath of we would like to read
   * \return the number of elements to the array
   */
  int get_length_of_table(const std::string& filename,
                          const std::string& datapath);

 public:
  // materialLibrary constructor
  MaterialLibrary();  //< empty constryctor

  /**
   * Constructor from file
   * \param filename path to file on disk, this file may be either in plaintext
   *                 or HDF5 format.
  */
  MaterialLibrary(const std::string& filename,
                  const std::string& datapath = "/materials");

  ~MaterialLibrary();  //< default destructor

  void renumber_mat();
  /**
   * \brief loads the pyne materials in map of name vs Material
    /// \param filename Path on disk to the HDF5 file.
    /// \param datapath Path to the materials in the file.
    /// \param protocol Flag for layout of material on disk.
  */
  void from_hdf5(const std::string& filename,
                 const std::string& datapath = "/materials",
                 const std::string& nucpath = "", int protocol = 1);
  /**
   * \brief loads the pyne materials in map of name vs Material
    /// \param filename Path on disk to the HDF5 file.
    /// \param datapath Path to the materials in the file.
    /// \param protocol Flag for layout of material on disk.
  */
  void from_hdf5(char* filename, char* datapath = "/materials",
                 char* nucpath = "", int protocol = 1);
  /**
   * \brief loads the pyne materials in map of name vs Material
    /// \param filename Path on disk to the json file.
  */
  void from_json(char* filename);
  /**
   * \brief loads the pyne materials in map of name vs Material
    /// \param filename Path on disk to the json file.
  */
  void from_json(const std::string& filename);
  void load_json(Json::Value json);
  Json::Value dump_json();
  void write_json(const std::string& filename);
  void write_json(char* filename);
  /**
   * Writes MaterialLibrary out to an HDF5 file.
   *  This happens according to protocol 1.
   *  \param filename Path on disk to the HDF5 file.
   *  \param datapath Path to the the material in the file.
   *  \param nucpath Path to the nuclides set in the file.
   *  \param row The index to read out, may be negative. Also note that this is
   * a
   *             float.  A value of -0.0 indicates that the material should be
   *             appended to the end of the dataset.
   *  \param chunksize The chunksize for all material data on disk.
  */
  void write_hdf5(char* filename, char* datapath = "/materials",
                  char* nucpath = "/nucid");
  /**
   * Writes MaterialLibrary out to an HDF5 file.
   *  This happens according to protocol 1.
   *  \param filename Path on disk to the HDF5 file.
   *  \param datapath Path to the the material in the file.
   *  \param nucpath Path to the nuclides set in the file.
   *  \param row The index to read out, may be negative. Also note that this is
   * a
   *             float.  A value of -0.0 indicates that the material should be
   *             appended to the end of the dataset.
   *  \param chunksize The chunksize for all material data on disk.
  */
  void write_hdf5(const std::string& filename,
                  const std::string& datapath = "/materials",
                  const std::string& nucpath = "/nucid");
  /**
   * \brief Merge a material library into the current one
   * \param mat_library pyne material library to merge
  */
  void merge(pyne::MaterialLibrary mat_lib);

  /**
   * \brief Add a material to the library
   * \param mat material to add
  */
  void add_material(pyne::Material mat);
  /**
   * \brief Add a material to the library
   * \param mat material to add
  */
  void add_material(pyne::Material mat, char* mat_name);
  /**
   * \brief Add a material to the library
   * \param mat material to add
  */
  void add_material(pyne::Material mat, const std::string& mat_name);
  /**
   * \brief replace a material to the library
   * \param num position of the material to replace
   * \param mat material to add
  */
  void replace(int num, pyne::Material mat);
  /**
   * \brief remove a material of the Library
   * \param mat material to remove
  */
  void del_material(pyne::Material mat);
  /**
   * \brief remove a material of the Library by name
   * \param mat_name name of the material to remove
  */
  void del_material(const std::string& mat_name);
  /**
   * \brief Get a material of the Library by name
   * \param mat_name name of the material to return
  */
  pyne::Material get_material(const std::string& mat_name);
  /**
   * \brief Get a material of the Library by name
   * \param num number of the material to return
  */
  pyne::Material get_material(int num);
  /**
   * \brief Get the material Library
   * \return std::map<std::string, pyne::MaterialLibrary>
  */
  pyne::mat_map get_mat_library() { return material_library; }
  /**
   * \brief Get the list of materials in the Library
   * \return std::set<std::string>
  */
  pyne::matname_set get_matlist() { return matlist; }
  /**
   * \brief Get the list of nuclides in the Library
   * \return std::set<int>
  */
  pyne::nuc_set get_nuclist() { return nuclist; }

  std::vector<std::string> name_order;

};  // end MaterialLibrary class header
}  // end of pyne namespace
#endif  // PYNE_MATERIAL_LIBRARY
