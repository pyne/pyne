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
