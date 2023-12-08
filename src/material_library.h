#ifndef PYNE_MATERIAL_LIBRARY
#define PYNE_MATERIAL_LIBRARY

#include <stdlib.h>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#if !defined(JSON_IS_AMALGAMATION)
#define JSON_IS_AMALGAMATION
#endif

#ifndef PYNE_IS_AMALGAMATED
#include "material.h"
#endif

namespace pyne {

typedef std::shared_ptr<pyne::Material> shr_mat_ptr;
typedef std::unordered_map<std::string, shr_mat_ptr> mat_map;
typedef std::set<std::string> matname_set;
typedef std::set<int> nuc_set;

class MaterialLibrary {
 public:
  // materialLibrary constructor
  MaterialLibrary();  //< empty constryctor

  /**
   * \brief Constructor from file
   * \param filename path to file on disk, this file may be either in
            plaintext or HDF5 format.
   * \param datapath Path to the materials in the
   * file.
   */
  MaterialLibrary(const std::string& fname, const std::string& datapath);

  ~MaterialLibrary(){};  //< default destructor

  /**
   * \brief loads the pyne materials in map of name vs Material
   * \param filename Path on disk to the HDF5 file.
   * \param datapath Path to the materials in the file.
   * \param protocol Flag for layout of material on disk.
   */
  void from_hdf5(const std::string& filename, const std::string& datapath);
  /**
   * \brief loads the pyne materials in map of name vs Material
    /// \param filename Path on disk to the json file.
  */
  void from_json(const std::string& filename);
  void load_json(Json::Value json);
  Json::Value dump_json();
  void write_json(const std::string& filename);
  //! Writes the materials in an OpenMC XML format
  void write_openmc(const std::string& filename) const;

  /**
   * \brief Writes MaterialLibrary out to an HDF5 file.
            This happens according to protocol 1.
            Writting in a file already containing
            This might cause issue to
            read material already present in the datapath,
            and for the new materials if written in the
            same datapath.
   *  \param filename Path on disk to the HDF5 file.
   *  \param datapath Path to the the material in the file.
   *  \param h5_overwrite bool to decide if existing path should be overwritten
   */
  void write_hdf5(const std::string& filename,
                  const std::string& datapath = "/materials",
                  bool h5_overwrite = false) const;
  /**
   * \brief Writes this nucpath to an HDF5 file.
   *  This happens according to protocol 1.
   * \param db HDF5 id for the open HDF5 file.
   * \param nucpath Path to the nuclides list in the file.
   */
  void write_hdf5_nucpath(hid_t db, std::string nucpath) const;
  /**
   * \brief Merge a material library into the current one
   * \param mat_library pyne material library to merge
   */
  void merge(const pyne::MaterialLibrary& mat_lib);
  /**
   * \brief Merge a material library into the current one
   * \param mat_library pointer to the pyne material library to merge
   */
  void merge(pyne::MaterialLibrary* mat_lib);
  /**
   * \brief Add a material to the library
   * \param mat material to add
   */
  void add_material(pyne::Material mat);
  /**
   * \brief Add a material to the library
   * \param mat_name name of the material to add(will overwrite material name
            if it already has one)
   * \param mat material to add
  */
  void add_material(const std::string& mat_name, const pyne::Material& mat);
  /**
   * \brief remove a material of the Library by name
   * \param mat_name name of the material to remove
   */
  void del_material(const std::string& mat_name);
  /**
   * \brief remove a material of the Library by name
   * \param mat_name name of the material to remove
   */
  pyne::Material get_material(const std::string& mat_name) const;
  /**
   * \brief Get a material of the Library by name
   * \param mat_name name of the material to return
   */
  pyne::shr_mat_ptr get_material_ptr(const std::string& mat_name) const;
  /**
 * \brief Return a material material number, ensure it exist if it does not
          exist build it using the material number, if material number is not
          defined, define it accordingly to the material library
 * \param mat material to check
 * \return std::string
 */
  std::string ensure_material_name_and_number(pyne::Material& mat) const;
  /**
 * \brief Return a material material number, ensure it exist as an int, if it
          exist as a string convert it into int and update the material
 * \param mat material to check
 * \return int
 */
  int ensure_material_number(pyne::Material& mat) const;
  /**
   * \brief Get the material library itself
   * \return std::map<std::string, pyne::Material>
   */
  pyne::mat_map get_mat_library() const { return material_library; }
  /**
   * \brief Get the list of materials in the Library
   * \return std::set<std::string>
   */
  pyne::matname_set get_keylist() const { return keylist; }
  /**
   * \brief Get the list of nuclides in the Library
   * \return std::set<int>
   */
  pyne::nuc_set get_nuclist() const { return nuclist; }

  /** Set of method allowing MaterialLibrary to behave like a undordered_map,
   * pointing directly to the nested material_library
   **/
  typedef mat_map::iterator iterator;
  typedef mat_map::const_iterator const_iterator;
  iterator begin() { return material_library.begin(); }
  iterator end() { return material_library.end(); }
  const_iterator begin() const { return material_library.begin(); }
  const_iterator end() const { return material_library.end(); }
  const_iterator cbegin() const { return material_library.cbegin(); }
  const_iterator cend() const { return material_library.cend(); }

  std::size_t size() const { return material_library.size(); }
  bool empty() const { return material_library.empty(); }
  std::size_t count(std::string mat_name) const {
    return material_library.count(mat_name);
  }
  std::shared_ptr<Material>& operator[](const std::string& k) {
    return material_library[k];
  }
  std::shared_ptr<Material>& operator[](std::string&& k) {
    return material_library[k];
  }

 private:
  /**
   * \brief Add the nuclide form the material to the list of nuclides
   * \param material from which add the nuclide to the list
   */
  void append_to_nuclist(const pyne::Material& mat);
  /**
   * \brief determines the length of an hdf5 data table
   * \param db HDF5 id for the open HDF5 file.
   * \param datapath we would like to read
   * \return the number of elements to the array
   */
  int get_length_of_table(hid_t db, const std::string& datapath) const;

  std::set<int>
      mat_number_set;   // list of material number in the library (ordered)
  matname_set keylist;  // list of masterial keys
  nuc_set nuclist;      // list of nuclides in the library
  // The actual library
  mat_map material_library;  // material library

};  // end MaterialLibrary class header
}  // namespace pyne
#endif  // PYNE_MATERIAL_LIBRARY
