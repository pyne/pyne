// Provides some HDF5 helper functionality in its own namespace

#if !defined(_H5_WRAP_)
#define _H5_WRAP_

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

#include "extra_types.h"


namespace h5wrap
{
  // Exceptions
  class HDF5BoundsError: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Index of point is out of bounds.  Cannot handle in HDF5 file.";
    };
  };


  class FileNotHDF5: public std::exception
  {
  public:
    FileNotHDF5(){};
    ~FileNotHDF5() throw () {};

    FileNotHDF5(std::string fname)
    {
      filename = fname;
    };

    virtual const char* what() const throw()
    {
      std::string FNH5str ("Not a valid HDF5 file: ");
      if (!filename.empty())
        FNH5str += filename;

      return (const char *) FNH5str.c_str();
    };

  private:
    std::string filename;
  };


  class GroupNotFound: public std::exception
  {
  public:
    GroupNotFound(){};
    ~GroupNotFound() throw () {};

    GroupNotFound(std::string fname, std::string gname)
    {
      filename = fname;
    };

    virtual const char* what() const throw()
    {
      std::string msg ("the group ");
      msg += groupname;
      msg += " not found in the file ";
      msg += filename;
      return (const char *) msg.c_str();
    };

  private:
    std::string filename;
    std::string groupname;
  };


  class PathNotFound: public std::exception
  {
  public:
    PathNotFound(){};
    ~PathNotFound() throw () {};

    PathNotFound(std::string fname, std::string pname)
    {
      filename = fname;
      path = pname;
    };

    virtual const char* what() const throw()
    {
      std::string msg ("the path ");
      msg += path;
      msg += " was not found in the HDF5 file ";
      msg += filename;
      return (const char *) msg.c_str();
    };

  private:
    std::string filename;
    std::string path;
  };



  // Read-in Functions
  template <typename T>
  T get_array_index(hid_t dset, int n, hid_t dtype=H5T_NATIVE_DOUBLE)
  {
    hsize_t count  [1] = {1};
    hsize_t offset [1] = {n};

    hid_t dspace = H5Dget_space(dset);
    hsize_t npoints = H5Sget_simple_extent_npoints(dspace);

    //Handle negative indices
    if (n < 0)
        offset[0] = offset[0] + npoints;

    //If still out of range we have a problem
    if (offset[0] < 0 || npoints <= offset[0])
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
  };


  // Conversion functions
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
    return cpp_set;
  };



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
    T mem_arr [arr_dims[0]];
    H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mem_arr);

    // Load new values into the vector
    cpp_vec.assign(mem_arr, mem_arr+arr_dims[0]);

    H5Dclose(dset);
    return cpp_vec;
  };


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
  };


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
  };



  // Classes
  template <typename T>
  class HomogenousTypeTable
  {
  public:
    HomogenousTypeTable(){};
    ~HomogenousTypeTable(){};
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
      T * col_buf;

      data.clear();
      for(int n = 0; n < shape[1]; n++)
      {
        // Make a compound data type of just this column
        col_type = H5Tcreate(H5T_COMPOUND, sizeof(T));
        H5Tinsert(col_type, cols[n].c_str(), n*sizeof(T), dtype);

        // Read in this column
        col_buf = new T [shape[0]];
        H5Dread(h5_set, col_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, col_buf);

        // save this column as a vector in out data map
        data[cols[n]] = std::vector<T>(col_buf, col_buf+shape[0]);
      };
    };

    // Metadata attribute
    std::string path;
    int shape [2];
    std::vector<std::string> cols;
    std::map<std::string, std::vector<T> > data;

    //
    // operator overloads
    //
    // index by column name (string)
    std::vector<T> operator[] (std::string col_name)
    {
      return data[col_name];
    };

    // index by int
    std::map<std::string, T> operator[] (int m)
    {
      std::map<std::string, T> row = std::map<std::string, T>();

      for(int n = 0; n < shape[1]; n++)
        row[cols[n]] = data[cols[n]][m];

      return row;
    };
  };


  /********************************/
  /*** Support for complex data ***/
  /********************************/
  hid_t _get_PYTABLES_COMPLEX128()
  {
    hid_t ct = H5Tcreate(H5T_COMPOUND, sizeof(extra_types::complex_t));
    H5Tinsert(ct, "r", HOFFSET(extra_types::complex_t, re), H5T_NATIVE_DOUBLE);
    H5Tinsert(ct, "i", HOFFSET(extra_types::complex_t, im), H5T_NATIVE_DOUBLE);
    return ct;
  };

  hid_t PYTABLES_COMPLEX128 = _get_PYTABLES_COMPLEX128();


  /*** Helper functions ***/
  bool path_exists(hid_t h5file, std::string path)
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
  };


// End namespace h5wrap
};



#endif

