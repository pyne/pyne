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

#include "H5Cpp.h"

#include "pyne.h"


namespace h5wrap
{
  // Exceptions
  class HDF5BoundsError: public std::exception
  {
    virtual const char* what() const throw()
    {
      return "Index of point is out of bounds.  Cannot read in from HDF5 File.";
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



  // Read-in Functions
  template <typename T>
  T get_array_index(H5::DataSet * ds, int n, H5::DataType dt = H5::PredType::NATIVE_DOUBLE)
  {
    H5::DataSpace array_space = (*ds).getSpace();

    hsize_t count  [1] = {1};
    hsize_t offset [1] = {n};

    //Handle negative indices
    if (n < 0)
        offset[0] = offset[0] + array_space.getSimpleExtentNpoints();

    //If still out of range we have a problem
    if (offset[0] < 0 || array_space.getSimpleExtentNpoints() <= offset[0])
        throw HDF5BoundsError();

    array_space.selectHyperslab(H5S_SELECT_SET, count, offset);

    //Set memmory hyperspace
    hsize_t dimsm[1] = {1};
    H5::DataSpace memspace(1, dimsm);

    hsize_t count_out  [1] = {1};
    hsize_t offset_out [1] = {0};

    memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    T data_out [1];
    (*ds).read(data_out, dt, memspace, array_space);

    return data_out[0];
  };


  // Conversion functions
  template <typename T>
  std::set<T> h5_array_to_cpp_set(H5::H5File * h5_file, std::string data_path, H5::DataType dt = H5::PredType::NATIVE_DOUBLE)
  {
    // Init
    std::set<T> cpp_set = std::set<T>();
    hsize_t arr_len[1];

    H5::DataSet h5_arr = (*h5_file).openDataSet(data_path);

    // Initilize to dataspace, to find the indices we are looping over
    H5::DataSpace arr_space = h5_arr.getSpace();
    int arr_dim = arr_space.getSimpleExtentDims(arr_len, NULL);

    // Allocate memory buffer    
    T * mem_arr = new T [arr_len[0]];

    // Read in data from file to memory
    h5_arr.read(mem_arr, dt);

    // Load new values into the set
    cpp_set.insert(&mem_arr[0], &mem_arr[arr_len[0]]);

    h5_arr.close();
    return cpp_set;
  };



  template <typename T>
  std::vector<T> h5_array_to_cpp_vector_1d(H5::H5File * h5_file, std::string data_path, H5::DataType dt = H5::PredType::NATIVE_DOUBLE)
  {
    // Init
    hsize_t arr_dims [1];
    H5::DataSet h5_arr = (*h5_file).openDataSet(data_path);

    // Initilize to dataspace, to find the indices we are looping over
    H5::DataSpace arr_space = h5_arr.getSpace();
    int arr_ndims = arr_space.getSimpleExtentDims(arr_dims, NULL);

    // Allocate memory buffer    
    T mem_arr [arr_dims[0]];

    // Read in data from file to memory
    h5_arr.read(mem_arr, dt);

    // Initialize vector 
    std::vector<T> cpp_vec;

    // Load new values into the vector
    cpp_vec.assign(mem_arr, mem_arr+arr_dims[0]);

    h5_arr.close();
    return cpp_vec;
  };


  template <typename T>
  std::vector< std::vector<T> > h5_array_to_cpp_vector_2d(H5::H5File * h5_file, std::string data_path, H5::DataType dt = H5::PredType::NATIVE_DOUBLE)
  {
    // Init
    hsize_t arr_dims [2];
    H5::DataSet h5_arr = (*h5_file).openDataSet(data_path);

    // Initilize to dataspace, to find the indices we are looping over
    H5::DataSpace arr_space = h5_arr.getSpace();
    int arr_ndims = arr_space.getSimpleExtentDims(arr_dims, NULL);

    // Allocate memory buffer    
    T mem_arr [arr_dims[0] * arr_dims[1]];

    // Read in data from file to memory
    // Have to read in as 1D array to get HDF5 and new keyword
    // to play nice with each other
    h5_arr.read(mem_arr, dt);

    // Initialize vector of vectors
    std::vector< std::vector<T> > cpp_vec (arr_dims[0], std::vector<T>(arr_dims[1]));

    // Load new values into the vector of vectors, using some indexing tricks
    for(int i = 0; i < arr_dims[0]; i++)
    {
        cpp_vec[i].assign(mem_arr+(i*arr_dims[1]), mem_arr+((i+1)*arr_dims[1]));
    };

    h5_arr.close();
    return cpp_vec;
  };


  template <typename T>
  std::vector< std::vector< std::vector<T> > > h5_array_to_cpp_vector_3d(H5::H5File * h5_file, std::string data_path, H5::DataType dt = H5::PredType::NATIVE_DOUBLE)
  {
    // Init
    hsize_t arr_dims [3];
    H5::DataSet h5_arr = (*h5_file).openDataSet(data_path);

    // Initilize to dataspace, to find the indices we are looping over
    H5::DataSpace arr_space = h5_arr.getSpace();
    int arr_ndims = arr_space.getSimpleExtentDims(arr_dims, NULL);

    // Allocate memory buffer    
    T mem_arr [arr_dims[0] * arr_dims[1] * arr_dims[2]];

    // Read in data from file to memory
    // Have to read in as 1D array to get HDF5 and new keyword
    // to play nice with each other
    h5_arr.read(mem_arr, dt);

    // Initialize vector of vectors of vectors
    std::vector< std::vector< std::vector<T> > > cpp_vec (arr_dims[0], std::vector< std::vector<T> >(arr_dims[1], std::vector<T>(arr_dims[2])));

    // Load new values into the vector of vectors of vectors, using some indexing tricks
    for(int i = 0; i < arr_dims[0]; i++)
    {
        for(int j = 0; j < arr_dims[1]; j++)
        {
            cpp_vec[i][j].assign(mem_arr+((i*arr_dims[1]*arr_dims[2]) + (j*arr_dims[2])), mem_arr+((i*arr_dims[1]*arr_dims[2]) + ((j+1)*arr_dims[2])));
        };
    };

    h5_arr.close();
    return cpp_vec;
  }



  // Classes
  template <typename T>
  class HomogenousTypeTable
  {
  public:
    HomogenousTypeTable(){};
    ~HomogenousTypeTable(){};
    HomogenousTypeTable(H5::H5File * h5_file, std::string data_path, H5::DataType dt = H5::PredType::NATIVE_DOUBLE)
    {
      // Init 
      H5::DataSet h5_set = (*h5_file).openDataSet(data_path);
      H5::DataSpace h5_space = h5_set.getSpace();
      H5::CompType h5_type = H5::CompType(h5_set);

      // set path
      path = data_path;

      // set shape
      shape[0] = h5_space.getSimpleExtentNpoints();
      shape[1] = h5_type.getNmembers();

      // set cols
      std::string * cols_buf = new std::string [shape[1]];
      for(int n = 0; n < shape[1]; n++)
      {
        cols_buf[n] = h5_type.getMemberName(n);
      };
      cols.assign(cols_buf, cols_buf+shape[1]);

      // set data
      H5::CompType col_type;
      T * col_buf;

      data.clear();
      for(int n = 0; n < shape[1]; n++)
      {
        // Make a compound data type of just this column
        col_type = H5::CompType(sizeof(T));
        col_type.insertMember(cols[n], 0, dt);

        // allocate space to read in this column
        col_buf = new T [shape[0]];

        // Read in this column
        h5_set.read(col_buf, col_type);

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

  // End HomogenousTypeTable
  };


  /********************************/
  /*** Support for complex data ***/
  /********************************/
  typedef struct {
    double re;   /* real part */
    double im;   /* imaginary part */
  } complex_t;

  H5::CompType _get_COMPLEX()
  {
    H5::CompType ct(sizeof(complex_t));
    ct.insertMember("real", HOFFSET(complex_t, re), H5::PredType::NATIVE_DOUBLE);
    ct.insertMember("imag", HOFFSET(complex_t, im), H5::PredType::NATIVE_DOUBLE);
    return ct;
  };

  H5::CompType COMPLEX = _get_COMPLEX();


// End namespace h5wrap
};



#endif

