/// \file h5wrap.h
/// \author Anthony Scopatz (scopatz\@gmail.com)
///
/// \brief Provides some HDF5 helper functionality in its own namespace

#ifndef PYNE_MRNAFG5GNZDNPCRPX3UCBZ5MFE
#define PYNE_MRNAFG5GNZDNPCRPX3UCBZ5MFE

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <exception>

#include "hdf5.h"

#ifndef PYNE_IS_AMALGAMATED
#include "extra_types.h"
#endif

//! Wrapper for standard HDF5 operations
namespace h5wrap
{
  /// Custom exception for HDF5 indexing errors.
  class HDF5BoundsError: public std::exception
  {
    /// returns error message.
    virtual const char* what() const throw()
    {
      return "Index of point is out of bounds.  Cannot handle in HDF5 file.";
    };
  };


  /// Custom exception for when an existing file is not in a valid HDF5 format.
  class FileNotHDF5: public std::exception
  {
  public:

    /// default constructor
    FileNotHDF5(){};

    /// default destructor
    ~FileNotHDF5() throw () {};

    /// constructor with the filename
    FileNotHDF5(std::string fname)
    {
      filename = fname;
    };

    /// helpful error message that includes the filename
    virtual const char* what() const throw()
    {
      std::string FNH5str ("Not a valid HDF5 file: ");
      if (!filename.empty())
        FNH5str += filename;

      const char* FNH5str_rtn = FNH5str.c_str();
      return FNH5str_rtn;
    };

  private:
    std::string filename; ///< the file which is not in HDF5 format.
  };


  /// Custom exception for when a group cannot be found in an HDF5 file.
  class GroupNotFound: public std::exception
  {
  public:

    /// default constructor
    GroupNotFound(){};

    /// default destructor
    ~GroupNotFound() throw () {};

    /// constructor with the filename and the groupname
    GroupNotFound(std::string fname, std::string gname)
    {
      filename = fname;
    };

    /// helpful error message that includes the filename and the groupname
    virtual const char* what() const throw()
    {
      std::string msg ("the group ");
      msg += groupname;
      msg += " not found in the file ";
      msg += filename;
      const char* msg_rtn = msg.c_str();
      return msg_rtn;
    };

  private:
    std::string filename;   ///< the HDF5 file
    std::string groupname;  ///< the group in the hierarchy
  };

  /// Custom exception for when a path is not found in an HDF5 file
  class PathNotFound: public std::exception
  {
  public:

    /// default constructor
    PathNotFound(){};

    /// default destructor
    ~PathNotFound() throw () {};

    /// constructor with the filename and the pathname
    PathNotFound(std::string fname, std::string pname)
    {
      filename = fname;
      path = pname;
    };

    /// helpful error message that includes the filename and the pathname
    virtual const char* what() const throw()
    {
      std::string msg ("the path ");
      msg += path;
      msg += " was not found in the HDF5 file ";
      msg += filename;
      const char* msg_rtn = msg.c_str();
      return msg_rtn;
    };

  private:
    std::string filename; ///< the HDF5 file
    std::string path;     ///< the path in the file
  };



  // Read-in Functions

  /// Retrieves the \a nth index out of the dataset \a dset (which has an HDF5
  /// datatype \a dtype).  The value is returned as the C/C++ type given by \a T.
  template <typename T>
  T get_array_index(hid_t dset, int n, hid_t dtype=H5T_NATIVE_DOUBLE)
  {
    hsize_t count  [1] = {1};
    hsize_t offset [1] = {static_cast<hsize_t>(n)};

    hid_t dspace = H5Dget_space(dset);
    hsize_t npoints = H5Sget_simple_extent_npoints(dspace);

    //Handle negative indices
    if (n < 0)
        offset[0] = offset[0] + npoints;

    //If still out of range we have a problem
    if (npoints <= offset[0])
        throw HDF5BoundsError();

    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    //Set memmory hyperspace
    hsize_t dimsm[1] = {1};
    hid_t memspace = H5Screate_simple(1, dimsm, NULL);

    hsize_t count_out  [1] = {1};
    hsize_t offset_out [1] = {0};

    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL,
                                 count_out, NULL);

    T data_out [1];
    H5Dread(dset, dtype, memspace, dspace, H5P_DEFAULT, data_out);

    return data_out[0];
  }


  // Conversion functions

  /// Reads in data from an HDF5 file as a C++ set.  \a T should roughly match
  /// \a dtype.
  /// \param h5file HDF5 file id for an open file.
  /// \param data_path path to the data in the open file.
  /// \param dtype HDF5 data type for the data set at \a data_path.
  /// \return an in memory set of type \a T.
  template <typename T>
  std::set<T> h5_array_to_cpp_set(hid_t h5file, std::string data_path, hid_t dtype=H5T_NATIVE_DOUBLE)
  {
    std::set<T> cpp_set = std::set<T>();
    hsize_t arr_len[1];
    hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

    // Initilize to dataspace, to find the indices we are looping over
    hid_t arr_space = H5Dget_space(dset);
    int arr_dim = H5Sget_simple_extent_dims(arr_space, arr_len, NULL);

    // Read in data from file to memory
    T * mem_arr = new T [arr_len[0]];
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

    // Load new values into the set
    cpp_set.insert(&mem_arr[0], &mem_arr[arr_len[0]]);

    H5Dclose(dset);

    delete[] mem_arr;
    return cpp_set;
  }


  /// Reads in data from an HDF5 file as a 1 dimiensional vector.  \a T should roughly
  /// match \a dtype.
  /// \param h5file HDF5 file id for an open file.
  /// \param data_path path to the data in the open file.
  /// \param dtype HDF5 data type for the data set at \a data_path.
  /// \return an in memory 1D vector of type \a T.
  template <typename T>
  std::vector<T> h5_array_to_cpp_vector_1d(hid_t h5file, std::string data_path,
                                           hid_t dtype=H5T_NATIVE_DOUBLE)
  {
    std::vector<T> cpp_vec;
    hsize_t arr_dims [1];
    hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

    // Initilize to dataspace, to find the indices we are looping over
    hid_t arr_space = H5Dget_space(dset);
    int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

    // Read in data from file to memory
    T* mem_arr = new T [arr_dims[0]];
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

    // Load new values into the vector
    cpp_vec.assign(mem_arr, mem_arr+arr_dims[0]);

    H5Dclose(dset);
    return cpp_vec;
  }


  /// Reads in data from an HDF5 file as a 2 dimiensional vector.  \a T should roughly
  /// match \a dtype.
  /// \param h5file HDF5 file id for an open file.
  /// \param data_path path to the data in the open file.
  /// \param dtype HDF5 data type for the data set at \a data_path.
  /// \return an in memory 2D vector of type \a T.
  template <typename T>
  std::vector< std::vector<T> > h5_array_to_cpp_vector_2d(hid_t h5file, std::string data_path,
                                                          hid_t dtype=H5T_NATIVE_DOUBLE)
  {
    hsize_t arr_dims [2];
    hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

    // Initilize to dataspace, to find the indices we are looping over
    hid_t arr_space = H5Dget_space(dset);
    int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

    // Read in data from file to memory
    // Have to read in as 1D array to get HDF5 and new keyword
    // to play nice with each other
    T mem_arr [arr_dims[0] * arr_dims[1]];
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

    // Load new values into the vector of vectors, using some indexing tricks
    std::vector< std::vector<T> > cpp_vec (arr_dims[0], std::vector<T>(arr_dims[1]));
    for(int i = 0; i < arr_dims[0]; i++)
    {
        cpp_vec[i].assign(mem_arr+(i*arr_dims[1]), mem_arr+((i+1)*arr_dims[1]));
    };

    H5Dclose(dset);
    return cpp_vec;
  }


  /// Reads in data from an HDF5 file as a 3 dimiensional vector.  \a T should roughly
  /// match \a dtype.
  /// \param h5file HDF5 file id for an open file.
  /// \param data_path path to the data in the open file.
  /// \param dtype HDF5 data type for the data set at \a data_path.
  /// \return an in memory 3D vector of type \a T.
  template <typename T>
  std::vector< std::vector< std::vector<T> > > h5_array_to_cpp_vector_3d(hid_t h5file,
                                                  std::string data_path,
                                                  hid_t dtype=H5T_NATIVE_DOUBLE)
  {
    hsize_t arr_dims [3];
    hid_t dset = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);

    // Initilize to dataspace, to find the indices we are looping over
    hid_t arr_space = H5Dget_space(dset);
    int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

    // Read in data from file to memory
    // Have to read in as 1D array to get HDF5 and new keyword
    // to play nice with each other
    T mem_arr [arr_dims[0] * arr_dims[1] * arr_dims[2]];
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

    // Load new values into the vector of vectors of vectors, using some indexing tricks
    std::vector< std::vector< std::vector<T> > > cpp_vec (arr_dims[0], std::vector< std::vector<T> >(arr_dims[1], std::vector<T>(arr_dims[2])));
    for(int i = 0; i < arr_dims[0]; i++)
    {
        for(int j = 0; j < arr_dims[1]; j++)
        {
            cpp_vec[i][j].assign(mem_arr+((i*arr_dims[1]*arr_dims[2]) + (j*arr_dims[2])), mem_arr+((i*arr_dims[1]*arr_dims[2]) + ((j+1)*arr_dims[2])));
        };
    };

    H5Dclose(dset);
    return cpp_vec;
  }



  // Classes
  /// A class representing a high-level table contruct whose columns all have the same
  /// type \a T in C/C++ (and the analogous type in HDF5).
  template <typename T>
  class HomogenousTypeTable
  {
  public:

    /// default constructor
    HomogenousTypeTable(){};

    /// default destructor
    ~HomogenousTypeTable(){};

    /// Constructor to load in data upon initialization.  \a T should roughly
    /// match \a dtype.
    /// \param h5file HDF5 file id for an open file.
    /// \param data_path path to the data in the open file.
    /// \param dtype HDF5 data type for the data set at \a data_path.
    HomogenousTypeTable(hid_t h5file, std::string data_path, hid_t dtype=H5T_NATIVE_DOUBLE)
    {
      hid_t h5_set = H5Dopen2(h5file, data_path.c_str(), H5P_DEFAULT);
      hid_t h5_space = H5Dget_space(h5_set);
      hid_t h5_type = H5Dget_type(h5_set);

      // set path
      path = data_path;

      // set shape
      shape[0] = H5Sget_simple_extent_npoints(h5_space);
      shape[1] = H5Tget_nmembers(h5_type);

      // set cols
      std::string * cols_buf = new std::string [shape[1]];
      for(int n = 0; n < shape[1]; n++)
        cols_buf[n] = H5Tget_member_name(h5_type, n);
      cols.assign(cols_buf, cols_buf+shape[1]);

      // set data
      hid_t col_type;
      T * col_buf = new T [shape[0]];

      data.clear();
      for(int n = 0; n < shape[1]; n++)
      {
        // Make a compound data type of just this column
        col_type = H5Tcreate(H5T_COMPOUND, sizeof(T));
        H5Tinsert(col_type, cols[n].c_str(), 0, dtype);

        // Read in this column
        H5Dread(h5_set, col_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, col_buf);

        // save this column as a vector in out data map
        data[cols[n]] = std::vector<T>(col_buf, col_buf+shape[0]);
      };
      delete[] col_buf;
    };

    // Metadata attributes
    std::string path; ///< path in file to the data
    int shape [2];    ///< table shape, rows x columns.
    std::vector<std::string> cols;  ///< column names
    /// mapping from column names to column data
    std::map<std::string, std::vector<T> > data;

    //
    // operator overloads
    //
    /// index into the table by column name (string)
    std::vector<T> operator[] (std::string col_name)
    {
      return data[col_name];
    };

    /// index into the table by row
    std::map<std::string, T> operator[] (int m)
    {
      std::map<std::string, T> row = std::map<std::string, T>();

      for(int n = 0; n < shape[1]; n++)
        row[cols[n]] = data[cols[n]][m];

      return row;
    };
  };


  /// Create an HDF5 data type for complex 128 bit data, which happens to match the
  /// complex data type that is used by PyTables ^_~.
  inline hid_t _get_PYTABLES_COMPLEX128()
  {
    hid_t ct = H5Tcreate(H5T_COMPOUND, sizeof(xd_complex_t));
    H5Tinsert(ct, "r", HOFFSET(xd_complex_t, re), H5T_NATIVE_DOUBLE);
    H5Tinsert(ct, "i", HOFFSET(xd_complex_t, im), H5T_NATIVE_DOUBLE);
    return ct;
  }

  /// The HDF5 id for a complex data type compatible with PyTables generated data.
  static hid_t PYTABLES_COMPLEX128 = _get_PYTABLES_COMPLEX128();


  /// Determines if a path exists in an hdf5 file.
  /// \param h5file HDF5 file id for an open file.
  /// \param path path to the data in the open file.
  /// \return true or false
  inline bool path_exists(hid_t h5file, std::string path)
  {
    bool rtn = false;
    hid_t ds = H5Dopen2(h5file, path.c_str(), H5P_DEFAULT);
    if (0 <= ds)
    {
      rtn = true;
      H5Dclose(ds);
    }
    else
    {
      hid_t grp = H5Gopen2(h5file, path.c_str(), H5P_DEFAULT);
      if (0 <= grp)
      {
        rtn = true;
        H5Gclose(grp);
      }
    }
    return rtn;
  }


// End namespace h5wrap
}



#endif
