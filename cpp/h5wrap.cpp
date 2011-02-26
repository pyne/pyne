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
T h5wrap::h5_array_to_cpp_set(H5::DataSet * h5_arr, std::set<T> * cpp_set, H5::DataType dt)
{
    // Init indexes
    T set_element;
    hsize_t arr_len[1];

    // clear out values currently in the set
    (*cpp_set).clear();

    // Initilize to dataspace, to find the indices we are looping over
    H5::DataSpace arr_space = (*h5_arr).getSpace();
    int arr_dim = arr_space.getSimpleExtentDims(arr_len, NULL);

    // Iterate over the elements of the array, adding them to the set.
    for(int n = 0; n < arr_len[0]; n++)
    {
        set_element = h5wrap::get_array_index<T>(h5_arr, n, dt);
        (*cpp_set).insert(set_element);
    };

};

template int h5wrap::h5_array_to_cpp_set(H5::DataSet *, std::set<int> *, H5::DataType);
template double h5wrap::h5_array_to_cpp_set(H5::DataSet *, std::set<double> *, H5::DataType);
