// Provides some HDF5 helper functionality in its own namespace

#if !defined(_H5_WRAP_)
#define _H5_WRAP_

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <exception>

#include "H5Cpp.h"

#include "bright.h"

namespace h5wrap
{
    // Read-in Functions
    template <typename T>
    T get_array_index(H5::DataSet *, int, H5::DataType = H5::PredType::NATIVE_DOUBLE);


    // Conversion functions
    template <typename T>
    T h5_array_to_cpp_set(H5::DataSet *, std::set<T> *, H5::DataType = H5::PredType::NATIVE_DOUBLE);


    // Exceptions
    class HDF5BoundsError: public std::exception
    {
      virtual const char* what() const throw()
      {
        return "Index of point is out of bounds.  Cannot read in from HDF5 File.";
      };
    };

};



#endif
