#ifndef PYNE_MYEOEP7PSJFDTKDNQ4V3BZDUQM
#define PYNE_MYEOEP7PSJFDTKDNQ4V3BZDUQM
/// \file endf.h
///
/// \brief Pure c++ endf parser
///
/// This c++ endf parser stores a library in the endf_library class. Evaluation
/// is then done using the interpolation and related functions

#include <sstream>
#include <fstream>
#include <list>
#include <map>

#ifndef PYNE_IS_AMALGAMATED
  #include "pyne.h"
  #include "endf_mt.h"
#endif

namespace pyne
{
namespace endf
{
  /// struct used for indexing endf library
  typedef struct endf_id {
    int mat;
    int mf;
    int mt;
    friend bool operator <(const endf_id &lhs, const endf_id &rhs);
    friend bool operator ==(const endf_id &lhs, const endf_id &rhs);
  } endf_id;

  bool operator <(const endf_id &lhs, const endf_id &rhs);
  bool operator ==(const endf_id &lhs, const endf_id &rhs);
  /// utility function to convert an mt_base into an id
  endf_id make_endf_id(mt_base input);

  /// Class for storing raw endf data
  class library {
  public:
      ~library();

      std::map<endf_id , mt_base*> contents;///< library data

      /// Access to library data by id struct
      template<typename T> T get (endf_id comp);
      /// Access to library data by mat number, mf, and mt number
      template<typename T> T get (int mat, int mf, int mt);
      /// Access to library data by mf and mt number returns a vector of
      /// structs of the templated type
      template<typename T> std::vector<T> getl (int mf, int mt);

      std::vector<std::vector<int> > get_content_list();
      /// add data in the file to this library
      void read_endf(std::string filenm);
  };

}
}

#endif
