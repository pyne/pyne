# Cython imports
from libcpp.utility cimport pair
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.vector cimport vector as cpp_vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free

# Python Imports
#cimport collections
import collections

cimport numpy as np
import numpy as np

# Local imports
cimport std
cimport extra_types

#
# Map conversions
#

# <int, int> conversions

cdef cpp_map[int, int] dict_to_map_int_int(dict pydict):
    cdef cpp_map[int, int] cppmap = cpp_map[int, int]()
    for key, value in pydict.items():
        cppmap[key] = value
    return cppmap

cdef dict map_to_dict_int_int(cpp_map[int, int] cppmap):
    pydict = {}
    cdef cpp_map[int, int].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = deref(mapiter).second
        inc(mapiter)

    return pydict



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




# 1D Integer arrays

cdef cpp_vector[int] array_to_vector_1d_int(np.ndarray[np.int32_t, ndim=1] arr):
    cdef cpp_vector[int] vec = cpp_vector[int]()
    cdef Py_ssize_t n, N 

    # Get and reserve the size of the vector
    # prevents excessive resizing
    N = arr.shape[0]
    vec.reserve(N)

    # Loop through the array
    for n in range(N):
        vec.push_back(arr[n])

    return vec


cdef np.ndarray[np.int32_t, ndim=1] vector_to_array_1d_int(cpp_vector[int] vec):
    cdef np.ndarray[np.int32_t, ndim=1] arr
    cdef int n, N

    # Get and reserve the size of the array
    N = vec.size()
    arr = np.zeros((N,), dtype=np.int32) 

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





# 3D Float arrays

cdef cpp_vector[cpp_vector[cpp_vector[double]]] array_to_vector_3d_dbl(np.ndarray[np.float64_t, ndim=3] arr):
    cdef Py_ssize_t i, I, j, J, k, K

    # Get and reserve the size of the vector
    # prevents excessive resizing
    I = arr.shape[0]
    J = arr.shape[1]
    K = arr.shape[2]

    cdef cpp_vector[cpp_vector[cpp_vector[double]]] vec = cpp_vector[cpp_vector[cpp_vector[double]]](I, cpp_vector[cpp_vector[double]](J, cpp_vector[double](K)))

    # Loop through the array
    for i in range(I):
        for j in range(J):
            for k in range(K):
                vec[i][j][k] = arr[i][j][k]

    return vec


cdef np.ndarray[np.float64_t, ndim=3] vector_to_array_3d_dbl(cpp_vector[cpp_vector[cpp_vector[double]]] vec):
    cdef np.ndarray[np.float64_t, ndim=3] arr
    cdef int i, I, j, J, k, K

    # Get and reserve the size of the array
    I = vec.size()
    J = vec[0].size()
    K = vec[0][0].size()
    arr = np.zeros((I, J, K), dtype=np.float64) 

    # loop through the vector
    for i in range(I):
        for j in range(J):
            for k in range(K):
                arr[i][j][k] = vec[i][j][k]

    return arr





#
# Map-Vector Conversions
#

# {int: np.array()} 
cdef cpp_map[int, cpp_vector[double]] dict_to_map_int_array_to_vector_1d_dbl(dict pydict):
    cdef cpp_map[int, cpp_vector[double]] cppmap = cpp_map[int, cpp_vector[double]]()

    for key, value in pydict.items():
        cppmap[key] = array_to_vector_1d_dbl(value)

    return cppmap


cdef dict map_to_dict_int_vector_to_array_1d_dbl(cpp_map[int, cpp_vector[double]] cppmap):
    pydict = {}
    cdef cpp_map[int, cpp_vector[double]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = vector_to_array_1d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict




cdef cpp_map[int, cpp_vector[cpp_vector[double]]] dict_to_map_int_array_to_vector_2d_dbl(dict pydict):
    cdef cpp_map[int, cpp_vector[cpp_vector[double]]] cppmap = cpp_map[int, cpp_vector[cpp_vector[double]]]()

    for key, value in pydict.items():
        cppmap[key] = array_to_vector_2d_dbl(value)

    return cppmap


cdef dict map_to_dict_int_vector_to_array_2d_dbl(cpp_map[int, cpp_vector[cpp_vector[double]]] cppmap):
    pydict = {}
    cdef cpp_map[int, cpp_vector[cpp_vector[double]]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = vector_to_array_2d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict




cdef cpp_map[int, cpp_vector[cpp_vector[cpp_vector[double]]]] dict_to_map_int_array_to_vector_3d_dbl(dict pydict):
    cdef cpp_map[int, cpp_vector[cpp_vector[cpp_vector[double]]]] cppmap = cpp_map[int, cpp_vector[cpp_vector[cpp_vector[double]]]]()

    for key, value in pydict.items():
        cppmap[key] = array_to_vector_3d_dbl(value)

    return cppmap


cdef dict map_to_dict_int_vector_to_array_3d_dbl(cpp_map[int, cpp_vector[cpp_vector[cpp_vector[double]]]] cppmap):
    pydict = {}
    cdef cpp_map[int, cpp_vector[cpp_vector[cpp_vector[double]]]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = vector_to_array_3d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict



# {str: np.array()} 
cdef cpp_map[std.string, cpp_vector[double]] dict_to_map_str_array_to_vector_1d_dbl(dict pydict):
    cdef std.string s
    cdef cpp_map[std.string, cpp_vector[double]] cppmap = cpp_map[std.string, cpp_vector[double]]()

    for key, value in pydict.items():
        s = std.string(key)
        cppmap[s] = array_to_vector_1d_dbl(value)

    return cppmap


cdef dict map_to_dict_str_vector_to_array_1d_dbl(cpp_map[std.string, cpp_vector[double]] cppmap):
    pydict = {}
    cdef cpp_map[std.string, cpp_vector[double]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[(deref(mapiter).first).c_str()] = vector_to_array_1d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict



#
# Map-Map-Vector Conversions
#

# {int: {int: np.array()}}
cdef cpp_map[int, cpp_map[int, cpp_vector[double]]] dict_to_map_int_int_array_to_vector_1d_dbl(dict pydict):
    cdef cpp_map[int, cpp_map[int, cpp_vector[double]]] cppmap = cpp_map[int, cpp_map[int, cpp_vector[double]]]()

    for key, value in pydict.items():
        cppmap[key] = dict_to_map_int_array_to_vector_1d_dbl(value)

    return cppmap


cdef dict map_to_dict_int_int_vector_to_array_1d_dbl(cpp_map[int, cpp_map[int, cpp_vector[double]]] cppmap):
    pydict = {}
    cdef cpp_map[int, cpp_map[int, cpp_vector[double]]].iterator mapiter = cppmap.begin()

    while mapiter != cppmap.end():
        pydict[deref(mapiter).first] = map_to_dict_int_vector_to_array_1d_dbl(deref(mapiter).second)
        inc(mapiter)

    return pydict



#
# Proxy Classes
#

#
# --- Sets
#

# Int
cdef class SetIterInt(object):
    cdef void init(self, cpp_set[int] * set_ptr):
        cdef cpp_set[int].iterator * itn = <cpp_set[int].iterator *> malloc(sizeof(set_ptr.begin()))
        itn[0] = set_ptr.begin()
        self.iter_now = itn

        cdef cpp_set[int].iterator * ite = <cpp_set[int].iterator *> malloc(sizeof(set_ptr.end()))
        ite[0] = set_ptr.end()
        self.iter_end = ite
        
    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_set[int].iterator inow = deref(self.iter_now)
        cdef cpp_set[int].iterator iend = deref(self.iter_end)
        cdef int val = deref(inow)

        if inow != iend:
            pyval = val
        else:
            raise StopIteration    

        inc(deref(self.iter_now))
        return pyval


cdef class _SetInt:
    def __cinit__(self, new_set=True, bint free_set=True):
        # Decide how to init set, if at all
        if isinstance(new_set, _SetInt):
            self.set_ptr = (<_SetInt> new_set).set_ptr
        elif hasattr(new_set, '__iter__') or (hasattr(new_set, '__len__') and hasattr(new_set, '__getitem__')):
            self.set_ptr = new cpp_set[int]()
            for value in new_set:
                self.set_ptr.insert(value)
        elif bool(new_set):
            self.set_ptr = new cpp_set[int]()

        # Store free_set
        self._free_set = free_set
        
    def __dealloc__(self):
        if self._free_set:
            del self.set_ptr

    def __contains__(self, value):
        if not isinstance(value, int):
            return False

        if 0 < self.set_ptr.count(value):
            return True
        else:
            return False

    def __len__(self):
        return self.set_ptr.size()

    def __iter__(self):
        cdef SetIterInt si = SetIterInt()
        si.init(self.set_ptr)
        return si

    # Add mutable interface

    def add(self, int value):
        self.set_ptr.insert(value)
        return 

    def discard(self, int value):
        if value in self:
            self.set_ptr.erase(value)
        return 



class SetInt(_SetInt, collections.MutableSet):
    """Wrapper class for C++ standard library sets of type <int>.
    Provides set like interface on the Python level.

    Parameters
    ----------
    new_set : bool or dict-like
        Boolean on whether to make a new set or not, or set-like object
        with keys and values which are castable to the appropriate type.
    free_set : bool
        Flag for whether the pointer to the C++ set should be deallocated
        when the wrapper is dereferenced.
    """
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "set([" + ", ".join([str(i) for i in self]) + "])"


# Str
cdef class SetIterStr(object):
    cdef void init(self, cpp_set[std.string] * set_ptr):
        cdef cpp_set[std.string].iterator * itn = <cpp_set[std.string].iterator *> malloc(sizeof(set_ptr.begin()))
        itn[0] = set_ptr.begin()
        self.iter_now = itn

        cdef cpp_set[std.string].iterator * ite = <cpp_set[std.string].iterator *> malloc(sizeof(set_ptr.end()))
        ite[0] = set_ptr.end()
        self.iter_end = ite
        
    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_set[std.string].iterator inow = deref(self.iter_now)
        cdef cpp_set[std.string].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = str(deref(inow).c_str())
        else:
            raise StopIteration    

        inc(deref(self.iter_now))
        return pyval


cdef class _SetStr:
    def __cinit__(self, new_set=True, bint free_set=True):
        cdef std.string s

        # Decide how to init set, if at all
        if isinstance(new_set, _SetStr):
            self.set_ptr = (<_SetStr> new_set).set_ptr
        elif hasattr(new_set, '__iter__') or (hasattr(new_set, '__len__') and hasattr(new_set, '__getitem__')):
            self.set_ptr = new cpp_set[std.string]()
            for value in new_set:
                s = std.string(value)
                self.set_ptr.insert(s)
        elif bool(new_set):
            self.set_ptr = new cpp_set[std.string]()

        # Store free_set
        self._free_set = free_set
        
    def __dealloc__(self):
        if self._free_set:
            del self.set_ptr

    def __contains__(self, value):
        cdef std.string s
        if isinstance(value, basestring):
            s = std.string(value)
        else:
            return False

        if 0 < self.set_ptr.count(s):
            return True
        else:
            return False

    def __len__(self):
        return self.set_ptr.size()

    def __iter__(self):
        cdef SetIterStr si = SetIterStr()
        si.init(self.set_ptr)
        return si

    # Add mutable interface

    def add(self, char * value):
        cdef std.string s
        s = std.string(value)
        self.set_ptr.insert(s)
        return 

    def discard(self, char * value):
        cdef std.string s
        if value in self:
            s = std.string(value)
            self.set_ptr.erase(s)
        return 


class SetStr(_SetStr, collections.Set):
    """Wrapper class for C++ standard library sets of type <string>.
    Provides set like interface on the Python level.

    Parameters
    ----------
    new_set : bool or dict-like
        Boolean on whether to make a new set or not, or set-like object
        with keys and values which are castable to the appropriate type.
    free_set : bool
        Flag for whether the pointer to the C++ set should be deallocated
        when the wrapper is dereferenced.
    """
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "set([" + ", ".join([repr(i) for i in self]) + "])"



#
# --- Maps
#

# (Str, Int)
cdef class MapIterStrInt(object):
    cdef void init(self, cpp_map[std.string, int] * map_ptr):
        cdef cpp_map[std.string, int].iterator * itn = <cpp_map[std.string, int].iterator *> malloc(sizeof(map_ptr.begin()))
        itn[0] = map_ptr.begin()
        self.iter_now = itn

        cdef cpp_map[std.string, int].iterator * ite = <cpp_map[std.string, int].iterator *> malloc(sizeof(map_ptr.end()))
        ite[0] = map_ptr.end()
        self.iter_end = ite
        
    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_map[std.string, int].iterator inow = deref(self.iter_now)
        cdef cpp_map[std.string, int].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = str(deref(inow).first.c_str())
        else:
            raise StopIteration    

        inc(deref(self.iter_now))
        return pyval


cdef class _MapStrInt:
    def __cinit__(self, new_map=True, bint free_map=True):
        cdef std.string s
        cdef pair[std.string, int] item

        # Decide how to init map, if at all
        if isinstance(new_map, _MapStrInt):
            self.map_ptr = (<_MapStrInt> new_map).map_ptr
        elif hasattr(new_map, 'items'):
            self.map_ptr = new cpp_map[std.string, int]()
            for key, value in new_map.items():
                s = std.string(key)
                item = pair[std.string, int](s, value)
                self.map_ptr.insert(item)
        elif hasattr(new_map, '__len__'):
            self.map_ptr = new cpp_map[std.string, int]()
            for i in new_map:
                s = std.string(i[0])
                item = pair[std.string, int](s, i[1])
                self.map_ptr.insert(item)
        elif bool(new_map):
            self.map_ptr = new cpp_map[std.string, int]()

        # Store free_map
        self._free_map = free_map
        
    def __dealloc__(self):
        if self._free_map:
            del self.map_ptr

    def __contains__(self, key):
        cdef std.string s
        if isinstance(key, str):
            s = std.string(key)
        else:
            return False

        if 0 < self.map_ptr.count(s):
            return True
        else:
            return False

    def __len__(self):
        return self.map_ptr.size()

    def __iter__(self):
        cdef MapIterStrInt mi = MapIterStrInt()
        mi.init(self.map_ptr)
        return mi

    def __getitem__(self, key):
        cdef std.string s
        if isinstance(key, basestring):
            s = std.string(key)
        else:
            raise TypeError("Only string keys are valid.")

        if 0 < self.map_ptr.count(s):
            return deref(self.map_ptr)[s]
        else:
            raise KeyError

    def __setitem__(self, char * key, int value):
        cdef std.string s = std.string(key)
        cdef pair[std.string, int] item = pair[std.string, int](s, value)
        self.map_ptr.insert(item)
        
    def __delitem__(self, char * key):
        cdef std.string s
        if key in self:
            s = std.string(key)
            self.map_ptr.erase(s)


class MapStrInt(_MapStrInt, collections.MutableMapping):
    """Wrapper class for C++ standard library maps of type <string, int>.
    Provides dictionary like interface on the Python level.

    Parameters
    ----------
    new_map : bool or dict-like
        Boolean on whether to make a new map or not, or dict-like object
        with keys and values which are castable to the appropriate type.
    free_map : bool
        Flag for whether the pointer to the C++ map should be deallocated
        when the wrapper is dereferenced.
    """
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "{" + ", ".join(["{0}: {1}".format(repr(key), value) for key, value in self.items()]) + "}"




# (Int, Str)
cdef class MapIterIntStr(object):
    cdef void init(self, cpp_map[int, std.string] * map_ptr):
        cdef cpp_map[int, std.string].iterator * itn = <cpp_map[int, std.string].iterator *> malloc(sizeof(map_ptr.begin()))
        itn[0] = map_ptr.begin()
        self.iter_now = itn

        cdef cpp_map[int, std.string].iterator * ite = <cpp_map[int, std.string].iterator *> malloc(sizeof(map_ptr.end()))
        ite[0] = map_ptr.end()
        self.iter_end = ite
        
    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_map[int, std.string].iterator inow = deref(self.iter_now)
        cdef cpp_map[int, std.string].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = int(deref(inow).first)
        else:
            raise StopIteration    

        inc(deref(self.iter_now))
        return pyval


cdef class _MapIntStr:
    def __cinit__(self, new_map=True, bint free_map=True):
        cdef std.string s
        cdef pair[int, std.string] item

        # Decide how to init map, if at all
        if isinstance(new_map, _MapIntStr):
            self.map_ptr = (<_MapIntStr> new_map).map_ptr
        elif hasattr(new_map, 'items'):
            self.map_ptr = new cpp_map[int, std.string]()
            for key, value in new_map.items():
                s = std.string(value)
                item = pair[int, std.string](key, s)
                self.map_ptr.insert(item)
        elif hasattr(new_map, '__len__'):
            self.map_ptr = new cpp_map[int, std.string]()
            for i in new_map:
                s = std.string(i[1])
                item = pair[int, std.string](i[0], s)
                self.map_ptr.insert(item)
        elif bool(new_map):
            self.map_ptr = new cpp_map[int, std.string]()

        # Store free_map
        self._free_map = free_map
        
    def __dealloc__(self):
        if self._free_map:
            del self.map_ptr

    def __contains__(self, key):
        if not isinstance(key, int):
            return False

        if 0 < self.map_ptr.count(key):
            return True
        else:
            return False

    def __len__(self):
        return self.map_ptr.size()

    def __iter__(self):
        cdef MapIterIntStr mi = MapIterIntStr()
        mi.init(self.map_ptr)
        return mi

    def __getitem__(self, key):
        if not isinstance(key, int):
            raise TypeError("Only integer keys are valid.")

        if 0 < self.map_ptr.count(key):
            return str(deref(self.map_ptr)[key].c_str())
        else:
            raise KeyError

    def __setitem__(self, int key, char * value):
        cdef std.string s = std.string(value)
        cdef pair[int, std.string] item = pair[int, std.string](key, s)
        self.map_ptr.insert(item)
        
    def __delitem__(self, int key):
        if key in self:
            self.map_ptr.erase(key)


class MapIntStr(_MapIntStr, collections.MutableMapping):
    """Wrapper class for C++ standard library maps of type <int, string>.
    Provides dictionary like interface on the Python level.

    Parameters
    ----------
    new_map : bool or dict-like
        Boolean on whether to make a new map or not, or dict-like object
        with keys and values which are castable to the appropriate type.
    free_map : bool
        Flag for whether the pointer to the C++ map should be deallocated
        when the wrapper is dereferenced.
    """
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "{" + ", ".join(["{0}: {1}".format(key, repr(value)) for key, value in self.items()]) + "}"






# (Int, Double)
cdef class MapIterIntDouble(object):
    cdef void init(self, cpp_map[int, double] * map_ptr):
        cdef cpp_map[int, double].iterator * itn = <cpp_map[int, double].iterator *> malloc(sizeof(map_ptr.begin()))
        itn[0] = map_ptr.begin()
        self.iter_now = itn

        cdef cpp_map[int, double].iterator * ite = <cpp_map[int, double].iterator *> malloc(sizeof(map_ptr.end()))
        ite[0] = map_ptr.end()
        self.iter_end = ite
        
    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_map[int, double].iterator inow = deref(self.iter_now)
        cdef cpp_map[int, double].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = int(deref(inow).first)
        else:
            raise StopIteration    

        inc(deref(self.iter_now))
        return pyval


cdef class _MapIntDouble:
    def __cinit__(self, new_map=True, bint free_map=True):
        cdef pair[int, double] item

        # Decide how to init map, if at all
        if isinstance(new_map, _MapIntDouble):
            self.map_ptr = (<_MapIntDouble> new_map).map_ptr
        elif hasattr(new_map, 'items'):
            self.map_ptr = new cpp_map[int, double]()
            for key, value in new_map.items():
                item = pair[int, double](key, value)
                self.map_ptr.insert(item)
        elif hasattr(new_map, '__len__'):
            self.map_ptr = new cpp_map[int, double]()
            for i in new_map:
                item = pair[int, double](i[0], i[1])
                self.map_ptr.insert(item)
        elif bool(new_map):
            self.map_ptr = new cpp_map[int, double]()

        # Store free_map
        self._free_map = free_map

    def __dealloc__(self):
        if self._free_map:
            del self.map_ptr

    def __contains__(self, key):
        if not isinstance(key, int):
            return False

        if 0 < self.map_ptr.count(key):
            return True
        else:
            return False

    def __len__(self):
        return self.map_ptr.size()

    def __iter__(self):
        cdef MapIterIntDouble mi = MapIterIntDouble()
        mi.init(self.map_ptr)
        return mi

    def __getitem__(self, key):
        if not isinstance(key, int):
            raise TypeError("Only integer keys are valid.")

        if 0 < self.map_ptr.count(key):
            return float(deref(self.map_ptr)[key])
        else:
            raise KeyError

    def __setitem__(self, int key, double value):
        cdef pair[int, double] item = pair[int, double](key, value)
        self.map_ptr.insert(item)
        
    def __delitem__(self, int key):
        if key in self:
            self.map_ptr.erase(key)


class MapIntDouble(_MapIntDouble, collections.MutableMapping):
    """Wrapper class for C++ standard library maps of type <int, double>.
    Provides dictionary like interface on the Python level.

    Parameters
    ----------
    new_map : bool or dict-like
        Boolean on whether to make a new map or not, or dict-like object
        with keys and values which are castable to the appropriate type.
    free_map : bool
        Flag for whether the pointer to the C++ map should be deallocated
        when the wrapper is dereferenced.
    """

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "{" + ", ".join(["{0}: {1}".format(key, value) for key, value in self.items()]) + "}"



#
# Map (Int, Complex)
#
cdef class MapIterIntComplex(object):
    cdef void init(self, cpp_map[int, extra_types.complex_t] * map_ptr):
        cdef cpp_map[int, extra_types.complex_t].iterator * itn = <cpp_map[int, extra_types.complex_t].iterator *> malloc(sizeof(map_ptr.begin()))
        itn[0] = map_ptr.begin()
        self.iter_now = itn

        cdef cpp_map[int, extra_types.complex_t].iterator * ite = <cpp_map[int, extra_types.complex_t].iterator *> malloc(sizeof(map_ptr.end()))
        ite[0] = map_ptr.end()
        self.iter_end = ite
        
    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_map[int, extra_types.complex_t].iterator inow = deref(self.iter_now)
        cdef cpp_map[int, extra_types.complex_t].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = int(deref(inow).first)
        else:
            raise StopIteration    

        inc(deref(self.iter_now))
        return pyval


cdef class _MapIntComplex:        
    def __cinit__(self, new_map=True, bint free_map=True):
        cdef pair[int, extra_types.complex_t] item
        cdef extra_types.complex_t v

        # Decide how to init map, if at all
        if isinstance(new_map, _MapIntComplex):
            self.map_ptr = (<_MapIntComplex> new_map).map_ptr
        elif hasattr(new_map, 'items'):
            self.map_ptr = new cpp_map[int, extra_types.complex_t]()
            for key, value in new_map.items():
                v = extra_types.complex_t()
                v.re = value.real
                v.im = value.imag
                item = pair[int, extra_types.complex_t](key, v)
                self.map_ptr.insert(item)
        elif hasattr(new_map, '__len__'):
            self.map_ptr = new cpp_map[int, extra_types.complex_t]()
            for i in new_map:
                v = extra_types.complex_t()
                v.re = i[1].real
                v.im = i[1].imag
                item = pair[int, extra_types.complex_t](i[0], v)
                self.map_ptr.insert(item)
        elif bool(new_map):
            self.map_ptr = new cpp_map[int, extra_types.complex_t]()

        # Store free_map
        self._free_map = free_map
        
    def __dealloc__(self):
        if self._free_map:
            del self.map_ptr

    def __contains__(self, key):
        if not isinstance(key, int):
            return False

        if 0 < self.map_ptr.count(key):
            return True
        else:
            return False

    def __len__(self):
        return self.map_ptr.size()

    def __iter__(self):
        cdef MapIterIntComplex mi = MapIterIntComplex()
        mi.init(self.map_ptr)
        return mi

    def __getitem__(self, key):
        cdef extra_types.complex_t value

        if not isinstance(key, int):
            raise TypeError("Only integer keys are valid.")

        if 0 < self.map_ptr.count(key):
            value = deref(self.map_ptr)[key]
            return complex(float(value.re), float(value.im))
        else:
            raise KeyError

    def __setitem__(self, int key, complex value):
        cdef extra_types.complex_t v = extra_types.complex_t()
        v.re = value.real
        v.im = value.imag
        cdef pair[int, extra_types.complex_t] item = pair[int, extra_types.complex_t](key, v)
        self.map_ptr.insert(item)
        
    def __delitem__(self, int key):
        if key in self:
            self.map_ptr.erase(key)



class MapIntComplex(_MapIntComplex, collections.MutableMapping):
    """Wrapper class for C++ standard library maps of type <int, complex>.
    Provides dictionary like interface on the Python level.

    Parameters
    ----------
    new_map : bool or dict-like
        Boolean on whether to make a new map or not, or dict-like object
        with keys and values which are castable to the appropriate type.
    free_map : bool
        Flag for whether the pointer to the C++ map should be deallocated
        when the wrapper is dereferenced.
    """

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "{" + ", ".join(["{0}: {1}".format(key, value) for key, value in self.items()]) + "}"
