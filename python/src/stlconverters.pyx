# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.vector cimport vector as cpp_vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

cimport numpy as np
import numpy as np

# Local imports
cimport std


#
# Map conversions
#

# <int, double> conversions

cdef cpp_map[int, double] dict_to_map_int_dbl(dict pydict):
    cdef cpp_map[int, double] cppmap = cpp_map[int, double]()
    for key, value in pydict.items():
        cppmap[key] = value
    return cppmap

cdef dict map_to_dict_int_dbl(cpp_map[int, double] cppmap):
    pydict = {}
    cdef cpp_map[int, double].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = deref(mapiter).second
        inc(mapiter)

    return pydict


# <string, int> conversions

cdef cpp_map[std.string, int] dict_to_map_str_int(dict pydict):
    cdef cpp_map[std.string, int] cppmap = cpp_map[std.string, int]()
    for key, value in pydict.items():
        cppmap[std.string(key)] = value
    return cppmap

cdef dict map_to_dict_str_int(cpp_map[std.string, int] cppmap):
    pydict = {}
    cdef cpp_map[std.string, int].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first.c_str()] = deref(mapiter).second
        inc(mapiter)

    return pydict


# <int, string> conversions

cdef cpp_map[int, std.string] dict_to_map_int_str(dict pydict):
    cdef cpp_map[int, std.string] cppmap = cpp_map[int, std.string]()
    for key, value in pydict.items():
        cppmap[key] = std.string(value)
    return cppmap


cdef dict map_to_dict_int_str(cpp_map[int, std.string] cppmap):
    pydict = {}
    cdef cpp_map[int, std.string].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = deref(mapiter).second.c_str()
        inc(mapiter)

    return pydict


# <string, double> conversions

cdef cpp_map[std.string, double] dict_to_map_str_dbl(dict pydict):
    cdef cpp_map[std.string, double] cppmap = cpp_map[std.string, double]()
    for key, value in pydict.items():
        cppmap[std.string(key)] = value
    return cppmap

cdef dict map_to_dict_str_dbl(cpp_map[std.string, double] cppmap):
    pydict = {}
    cdef cpp_map[std.string, double].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first.c_str()] = deref(mapiter).second
        inc(mapiter)

    return pydict






#
# Set conversions
#

# Integer sets

cdef cpp_set[int] py_to_cpp_set_int(set pyset):
    cdef cpp_set[int] cppset = cpp_set[int]()
    for item in pyset:
        cppset.insert(item)
    return cppset

cdef set cpp_to_py_set_int(cpp_set[int] cppset):
    pyset = set()
    cdef cpp_set[int].iterator setiter = cppset.begin()
    while setiter != cppset.end():
        pyset.add(deref(setiter))
        inc(setiter)
    return pyset


# String sets

cdef cpp_set[std.string] py_to_cpp_set_str(set pyset):
    cdef std.string s
    cdef cpp_set[std.string] cppset = cpp_set[std.string]()

    for item in pyset:
        s = std.string(item)
        cppset.insert(s)

    return cppset

cdef set cpp_to_py_set_str(cpp_set[std.string] cppset):
    pyset = set()
    cdef cpp_set[std.string].iterator setiter = cppset.begin()

    while setiter != cppset.end():
        pyset.add(deref(setiter).c_str())
        inc(setiter)

    return pyset



#
# Vector conversions
#

# 1D Float arrays

cdef cpp_vector[double] array_to_vector_1d_dbl(np.ndarray[np.float64_t, ndim=1] arr):
    cdef cpp_vector[double] vec = cpp_vector[double]()
    cdef Py_ssize_t n, N 

    # Get and reserve the size of the vector
    # prevents excessive resizing
    N = arr.shape[0]
    vec.reserve(N)

    # Loop through the array
    for n in range(N):
        vec.push_back(arr[n])

    return vec


cdef np.ndarray[np.float64_t, ndim=1] vector_to_array_1d_dbl(cpp_vector[double] vec):
    cdef np.ndarray[np.float64_t, ndim=1] arr
    cdef int n, N

    # Get and reserve the size of the array
    N = vec.size()
    arr = np.zeros((N,), dtype=np.float64) 

    # loop through the vector
    for n in range(N):
        arr[n] = vec[n]

    return arr




# 2D Float arrays

cdef cpp_vector[cpp_vector[double]] array_to_vector_2d_dbl(np.ndarray[np.float64_t, ndim=2] arr):
    cdef Py_ssize_t i, I, j, J 

    # Get and reserve the size of the vector
    # prevents excessive resizing
    I = arr.shape[0]
    J = arr.shape[1]

    cdef cpp_vector[cpp_vector[double]] vec = cpp_vector[cpp_vector[double]](I, cpp_vector[double](J))

    # Loop through the array
    for i in range(I):
        for j in range(J):
            vec[i][j] = arr[i][j]

    return vec


cdef np.ndarray[np.float64_t, ndim=2] vector_to_array_2d_dbl(cpp_vector[cpp_vector[double]] vec):
    cdef np.ndarray[np.float64_t, ndim=2] arr
    cdef int i, I, j, J

    # Get and reserve the size of the array
    I = vec.size()
    J = vec[0].size()
    arr = np.zeros((I, J), dtype=np.float64) 

    # loop through the vector
    for i in range(I):
        for j in range(J):
            arr[i][j] = vec[i][j]

    return arr





#
# Map-Vector Conversions
#

# {int: np.array()} 
cdef cpp_map[int, cpp_vector[double]] dict_to_map_int_vector_to_array_1d_dbl(dict pydict):
    cdef cpp_map[int, cpp_vector[double]] cppmap = cpp_map[int, cpp_vector[double]]()

    for key, value in pydict.items():
        cppmap[key] = array_to_vector_1d_dbl(value)

    return cppmap


cdef dict map_to_dict_int_array_to_vector_1d_dbl(cpp_map[int, cpp_vector[double]] cppmap):
    pydict = {}
    cdef cpp_map[int, cpp_vector[double]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = vector_to_array_1d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict




cdef cpp_map[int, cpp_vector[cpp_vector[double]]] dict_to_map_int_vector_to_array_2d_dbl(dict pydict):
    cdef cpp_map[int, cpp_vector[cpp_vector[double]]] cppmap = cpp_map[int, cpp_vector[cpp_vector[double]]]()

    for key, value in pydict.items():
        cppmap[key] = array_to_vector_2d_dbl(value)

    return cppmap


cdef dict map_to_dict_int_array_to_vector_2d_dbl(cpp_map[int, cpp_vector[cpp_vector[double]]] cppmap):
    pydict = {}
    cdef cpp_map[int, cpp_vector[cpp_vector[double]]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = vector_to_array_2d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict




#
# Map-Vector Conversions
#

# {int: {int: np.array()}}
cdef cpp_map[int, cpp_map[int, cpp_vector[double]]] dict_to_map_int_int_vector_to_array_1d_dbl(dict pydict):
    cdef cpp_map[int, cpp_map[int, cpp_vector[double]]] cppmap = cpp_map[int, cpp_map[int, cpp_vector[double]]]()

    for key, value in pydict.items():
        cppmap[key] = dict_to_map_int_vector_to_array_1d_dbl(value)

    return cppmap


cdef dict map_to_dict_int_int_array_to_vector_1d_dbl(cpp_map[int, cpp_map[int, cpp_vector[double]]] cppmap):
    pydict = {}
    cdef cpp_map[int, cpp_map[int, cpp_vector[double]]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = map_to_dict_int_array_to_vector_1d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict
