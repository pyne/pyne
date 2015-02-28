#ifndef PYNE_FYAUH7LGDRBBHMVWZJQBKGYECI
#define PYNE_FYAUH7LGDRBBHMVWZJQBKGYECI

#include <fstream>
#include <list>
#include <vector>

namespace pyne
{
  namespace endf {
    /*****************************************************************************/
    /******************* Structs for basic ENDF datatypes ************************/
    /*****************************************************************************/


    /// a struct matching the ENDF control structure.
    typedef struct control {
      double c1;  ///< First value in control struct ZZZAAA.M for HEAD [double]
      double c2;  ///< Second value in control struct AWR for HEAD [double]
      int l1; ///< first int in control struct
      int l2; ///< second int in control struct
      int n1; ///< third int in control struct
      int n2; ///< fourth int in control struct
      int mat; ///< material id
      int mf; ///< material file number
      int mt; ///< reaction number
    } control;

    /// a struct mathing the ENDF LIST structure
    typedef struct list {
      double c1;  ///< First value in control struct ZZZAAA.M for HEAD [double]
      double c2;  ///< Second value in control struct AWR for HEAD [double]
      int l1; ///< first int in control struct
      int l2; ///< second int in control struct
      int npl; ///< number of data points
      int n2; ///< fourth int in control struct
      int mat; ///< material id
      int mf; ///< material file number
      int mt; ///< reaction number
      std::vector<double> data; ///< data as doubles
    } list;


    typedef struct tab1 {
      double c1;  ///< First value in control struct ZZZAAA.M for HEAD [double]
      double c2;  ///< Second value in control struct AWR for HEAD [double]
      int l1; ///< first int in control struct
      int l2; ///< second int in control struct
      int nr; ///< number of data points
      int np; ///< fourth int in control struct
      int mat; ///< material id
      int mf; ///< material file number
      int mt; ///< reaction number
      std::vector<int> nbt; ///<
      std::vector<int> intn; ///< vector of interpolation type used
      std::vector<double> x;
      std::vector<double> y;
    } tab1;

    typedef struct tab2 {
      double c1;  ///< [double]
      double c2;  ///< [double]
      int l1; ///< first int in control struct
      int l2; ///< second int in control struct
      int nr; ///< number of data points
      int nz; ///< fourth int in control struct
      int mat; ///< material id
      int mf; ///< material file number
      int mt; ///< reaction number
      std::vector<int> nbt; ///<
      std::vector<int> intn; ///< vector of interpolation table
    } tab2;

    typedef struct intg {
      int ndigit;
      std::vector<std::vector<int> > intmat; ///< a multidimensional array maybe like this
    } intg;

    /// Basic ENDF record readers
    list read_list(std::ifstream &infile);
    control read_cont(std::ifstream &infile);
    tab1 read_tab1(std::ifstream &infile);
    tab2 read_tab2(std::ifstream &infile);
  }
}

#endif
