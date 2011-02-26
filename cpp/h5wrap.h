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

#include "bright.h"

namespace h5wrap
{
    // Read-in Functions
    template <typename T>
    T get_array_index(H5::DataSet *, int, H5::DataType = H5::PredType::NATIVE_DOUBLE);

    // Conversion functions
    template <typename T>
    std::set<T> h5_array_to_cpp_set(H5::H5File *, std::string, H5::DataType = H5::PredType::NATIVE_DOUBLE);

    // Classes
    template <typename T>
    class HomogenousTypeTable
    {
    public:
        HomogenousTypeTable();
        HomogenousTypeTable(H5::H5File *, std::string, H5::DataType = H5::PredType::NATIVE_DOUBLE);

        // Metadata attribute
        std::string path;
        int shape [2];
        std::vector<std::string> cols;
        std::map<std::string, std::vector<T> > data;

        // operator overloads
        std::vector<T> operator[] (std::string);
        std::map<std::string, T> operator[] (int);
    };


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
