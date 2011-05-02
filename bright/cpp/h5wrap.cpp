// h5wrap.cpp

#include "h5wrap.h"


/* 
 *  Index into an hdf5 array.
 */

template <typename T>
T h5wrap::get_array_index(H5::DataSet * ds, int n, H5::DataType dt)
{
    H5::DataSpace array_space = (*ds).getSpace();

    hsize_t count  [1] = {1};    
    hsize_t offset [1] = {n};   

    //Handle negative indices
    if (n < 0)
        offset[0] = offset[0] + array_space.getSimpleExtentNpoints();

    //If still out of range we have a problem
    if (offset[0] < 0 || array_space.getSimpleExtentNpoints() <= offset[0])
        throw h5wrap::HDF5BoundsError();

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

template int h5wrap::get_array_index(H5::DataSet *, int, H5::DataType);
template double h5wrap::get_array_index(H5::DataSet *, int, H5::DataType);




/* 
 *  Convert hdf5 array to C++ set.
 */

template <typename T>
std::set<T> h5wrap::h5_array_to_cpp_set(H5::H5File * h5_file, std::string data_path, H5::DataType dt)
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

    // Close out data set
    h5_arr.close();

    // Return set
    return cpp_set;
};


template std::set<int> h5wrap::h5_array_to_cpp_set(H5::H5File *, std::string, H5::DataType);
template std::set<double> h5wrap::h5_array_to_cpp_set(H5::H5File *, std::string, H5::DataType);






/* 
 *  Convert hdf5 array to C++ vector.
 */

template <typename T>
std::vector<T> h5wrap::h5_array_to_cpp_vector_1d(H5::H5File * h5_file, std::string data_path, H5::DataType dt)
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

    // Close out data set
    h5_arr.close();

    // Return vector of vectors
    return cpp_vec;
};

template std::vector<int> h5wrap::h5_array_to_cpp_vector_1d(H5::H5File *, std::string, H5::DataType);
template std::vector<double> h5wrap::h5_array_to_cpp_vector_1d(H5::H5File *, std::string, H5::DataType);




/* 
 *  Convert hdf5 array to C++ vector of vectors.
 */

template <typename T>
std::vector< std::vector<T> > h5wrap::h5_array_to_cpp_vector_2d(H5::H5File * h5_file, std::string data_path, H5::DataType dt)
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

    // Close out data set
    h5_arr.close();

    // Return vector of vectors
    return cpp_vec;
};

template std::vector< std::vector<int> > h5wrap::h5_array_to_cpp_vector_2d(H5::H5File *, std::string, H5::DataType);
template std::vector< std::vector<double> > h5wrap::h5_array_to_cpp_vector_2d(H5::H5File *, std::string, H5::DataType);




/* 
 *  Convert hdf5 array to C++ vector of vectors of vectors.
 */

template <typename T>
std::vector< std::vector< std::vector<T> > > h5wrap::h5_array_to_cpp_vector_3d(H5::H5File * h5_file, std::string data_path, H5::DataType dt)
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

    // Close out data set
    h5_arr.close();

    // Return vector of vectors
    return cpp_vec;
};

template std::vector< std::vector< std::vector<int> > > h5wrap::h5_array_to_cpp_vector_3d(H5::H5File * h5_file, std::string data_path, H5::DataType dt);
template std::vector< std::vector< std::vector<double> > > h5wrap::h5_array_to_cpp_vector_3d(H5::H5File * h5_file, std::string data_path, H5::DataType dt);






/*
 * Classes
 */



/*
 * HomogenousTypeTable
 */

template <typename T>
h5wrap::HomogenousTypeTable<T>::HomogenousTypeTable()
{
};


template <typename T>
h5wrap::HomogenousTypeTable<T>::HomogenousTypeTable(H5::H5File * h5_file, std::string data_path, H5::DataType dt)
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


template <typename T>
std::vector<T> h5wrap::HomogenousTypeTable<T>::operator[](std::string col_name)
{
    return data[col_name];
};


template <typename T>
std::map<std::string, T> h5wrap::HomogenousTypeTable<T>::operator[](int m)
{
    // init row
    std::map<std::string, T> row = std::map<std::string, T>(); 

    // fill row values
    for(int n = 0; n < shape[1]; n++)
    {
        row[cols[n]] = data[cols[n]][m];
    };

    // return row map
    return row;
};

template class h5wrap::HomogenousTypeTable<int>;
template class h5wrap::HomogenousTypeTable<double>;
