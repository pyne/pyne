"""Generates cython wrapper classes and converter functions for standard library
containters to the associated python types.
"""


ctypes = {
    'str': 'std_string',
    'int': 'int',
    'uint': 'unsigned int',
    'float': 'float',
    'double': 'double',
    'complex': 'extra_types.complex_t',
    }

cytypes = {
    'str': 'char *',
    'int': 'int',
    'uint': 'int',
    'float': 'float',
    'double': 'float',
    'complex': 'object',
    }

pytypes = {
    'str': 'basestring',
    'int': 'int',
    'uint': 'int',
    'float': 'float',
    'double': 'float',
    'complex': 'complex',
    }

class_names = {
    'str': 'Str',
    'int': 'Int',
    'uint': 'UInt',
    'float': 'Float',
    'double': 'Double',
    'complex': 'Complex',
    }

func_names = {
    'str': 'str',
    'int': 'int',
    'uint': 'uint',
    'float': 'flt',
    'double': 'dbl',
    'complex': 'complex',
    }

human_names = {
    'str': 'string',
    'int': 'integer',
    'uint': 'unsigned integer',
    'float': 'float',
    'double': 'double',
    'complex': 'complex',
    }

c2py_exprs = {
    'str': 'str(<char *> {var}.c_str())',
    'int': '{var}',
    'uint': '<int> {var}',
    'float': '{var}',
    'double': '{var}',
    'complex': 'complex(float({var}.re), float({var}.im))',
    }

py2c_exprs = {
    'str': 'std_string(<char *> {var})',
    'int': '{var}',
    'uint': '<unsigned int> {var}',
    'float': '<float> {var}',
    'double': '<double> {var}',
    'complex': 'py2c_complex({var})',
    }

testvals = {
    'str': ["Aha", "Take", "Me", "On"], 
    'int': [1, 42, -65, 18], 
    'uint': [1, 42, 65, 18],
    'float': [1.0, 42.42, -65.5555, 18],
    'double': [1.0, 42.42, -65.5555, 18],
    'complex': [1.0, 42+42j, -65.55-1j, 0.18j],
    }

#
# Sets
#

_pyxset = '''# Set{clsname}
cdef class SetIter{clsname}(object):
    cdef void init(self, cpp_set[{ctype}] * set_ptr):
        cdef cpp_set[{ctype}].iterator * itn = <cpp_set[{ctype}].iterator *> malloc(sizeof(set_ptr.begin()))
        itn[0] = set_ptr.begin()
        self.iter_now = itn

        cdef cpp_set[{ctype}].iterator * ite = <cpp_set[{ctype}].iterator *> malloc(sizeof(set_ptr.end()))
        ite[0] = set_ptr.end()
        self.iter_end = ite

    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_set[{ctype}].iterator inow = deref(self.iter_now)
        cdef cpp_set[{ctype}].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = {iterval}
        else:
            raise StopIteration

        inc(deref(self.iter_now))
        return pyval


cdef class _Set{clsname}:
    def __cinit__(self, new_set=True, bint free_set=True):
        cdef {ctype} s

        # Decide how to init set, if at all
        if isinstance(new_set, _Set{clsname}):
            self.set_ptr = (<_Set{clsname}> new_set).set_ptr
        elif hasattr(new_set, '__iter__') or \
                (hasattr(new_set, '__len__') and
                hasattr(new_set, '__getitem__')):
            self.set_ptr = new cpp_set[{ctype}]()
            for value in new_set:
                s = {initval}
                self.set_ptr.insert(s)
        elif bool(new_set):
            self.set_ptr = new cpp_set[{ctype}]()

        # Store free_set
        self._free_set = free_set

    def __dealloc__(self):
        if self._free_set:
            del self.set_ptr

    def __contains__(self, value):
        cdef {ctype} s
        if isinstance(value, {pytype}):
            s = {initval}
        else:
            return False

        if 0 < self.set_ptr.count(s):
            return True
        else:
            return False

    def __len__(self):
        return self.set_ptr.size()

    def __iter__(self):
        cdef SetIter{clsname} si = SetIter{clsname}()
        si.init(self.set_ptr)
        return si

    def add(self, {cytype} value):
        cdef {ctype} v
        v = {initval}
        self.set_ptr.insert(v)
        return

    def discard(self, value):
        cdef {ctype} v
        if value in self:
            v = {initval}
            self.set_ptr.erase(v)
        return


class Set{clsname}(_Set{clsname}, collections.Set):
    """Wrapper class for C++ standard library sets of type <{humname}>.
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

'''
def genpyx_set(t):
    """Returns the pyx snippet for a set of type t."""
    iterval = c2py_exprs[t].format(var="deref(inow)")
    initval = py2c_exprs[t].format(var="value")
    return _pyxset.format(clsname=class_names[t], humname=human_names[t], 
                          ctype=ctypes[t], pytype=pytypes[t], cytype=cytypes[t], 
                          iterval=iterval, initval=initval)

_pxdset = """# Set{clsname}
cdef class SetIter{clsname}(object):
    cdef cpp_set[{ctype}].iterator * iter_now
    cdef cpp_set[{ctype}].iterator * iter_end
    cdef void init(SetIter{clsname}, cpp_set[{ctype}] *)

cdef class _Set{clsname}:
    cdef cpp_set[{ctype}] * set_ptr
    cdef public bint _free_set


"""
def genpxd_set(t):
    """Returns the pxd snippet for a set of type t."""
    return _pyxset.format(clsname=class_names[t], ctype=ctypes[t])


_testset = """# Set{clsname}
def test_set_{fncname}():
    s = conv.Set{clsname}()
    s.add({0})
    assert_true({0} in s)
    assert_true({2} not in s)

    s = conv.Set{clsname}([{0}, {1}, {2}])
    assert_true({1} in s)
    assert_true({3} not in s)

"""
def gentest_set(t):
    """Returns the test snippet for a set of type t."""
    return _testset.format(*[repr(i) for i in testvals[t]], clsname=class_names[t],
                           fncname=func_names[t])

#
# Maps
#
_pyxmap = '''# Map({tclsname}, {uclsname})
cdef class MapIter{tclsname}{uclsname}(object):
    cdef void init(self, cpp_map[{tctype}, {uctype}] * map_ptr):
        cdef cpp_map[{tctype}, {uctype}].iterator * itn = <cpp_map[{tctype}, {uctype}].iterator *> malloc(sizeof(map_ptr.begin()))
        itn[0] = map_ptr.begin()
        self.iter_now = itn

        cdef cpp_map[{tctype}, {uctype}].iterator * ite = <cpp_map[{tctype}, {uctype}].iterator *> malloc(sizeof(map_ptr.end()))
        ite[0] = map_ptr.end()
        self.iter_end = ite

    def __dealloc__(self):
        free(self.iter_now)
        free(self.iter_end)

    def __iter__(self):
        return self

    def __next__(self):
        cdef cpp_map[{tctype}, {uctype}].iterator inow = deref(self.iter_now)
        cdef cpp_map[{tctype}, {uctype}].iterator iend = deref(self.iter_end)

        if inow != iend:
            pyval = {iterkey}
        else:
            raise StopIteration

        inc(deref(self.iter_now))
        return pyval

cdef class _Map{tclsname}{uclsname}:
    def __cinit__(self, new_map=True, bint free_map=True):
        cdef pair[{tctype}, {uctype}] item

        # Decide how to init map, if at all
        if isinstance(new_map, _Map{tclsname}{uclsname}):
            self.map_ptr = (<_Map{tclsname}{uclsname}> new_map).map_ptr
        elif hasattr(new_map, 'items'):
            self.map_ptr = new cpp_map[{tctype}, {uctype}]()
            for key, value in new_map.items():
                item = pair[{tctype}, {uctype}]({initkey}, {initval})
                self.map_ptr.insert(item)
        elif hasattr(new_map, '__len__'):
            self.map_ptr = new cpp_map[{tctype}, {uctype}]()
            for key, value in new_map:
                item = pair[{tctype}, {uctype}]({initkey}, {initval})
                self.map_ptr.insert(item)
        elif bool(new_map):
            self.map_ptr = new cpp_map[{tctype}, {uctype}]()

        # Store free_map
        self._free_map = free_map

    def __dealloc__(self):
        if self._free_map:
            del self.map_ptr

    def __contains__(self, key):
        cdef {tctype} k
        if not isinstance(key, {tpytype}):
            return False
        k = {initkey}

        if 0 < self.map_ptr.count(k):
            return True
        else:
            return False

    def __len__(self):
        return self.map_ptr.size()

    def __iter__(self):
        cdef MapIter{tclsname}{uclsname} mi = MapIter{tclsname}{uclsname}()
        mi.init(self.map_ptr)
        return mi

    def __getitem__(self, key):
        cdef {tctype} k
        cdef {uctype} v

        if not isinstance(key, {tpytype}):
            raise TypeError("Only {thumname} keys are valid.")
        k = {initkey}

        if 0 < self.map_ptr.count(k):
            v = deref(self.map_ptr)[k]
            return {convval}
        else:
            raise KeyError

    def __setitem__(self, key, value):
        cdef pair[{tctype}, {uctype}] item = pair[{tctype}, {uctype}]({initkey}, {initval})
        self.map_ptr.insert(item)

    def __delitem__(self, key):
        cdef {tctype} k
        if key in self:
            k = {initkey}
            self.map_ptr.erase(k)


class Map{tclsname}{uclsname}(_Map{tclsname}{uclsname}, collections.MutableMapping):
    """Wrapper class for C++ standard library maps of type <{thumname}, {uhumname}>.
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
        return "{{" + ", ".join(["{{0}}: {{1}}".format(key, value) for key, value in self.items()]) + "}}"

'''
def genpyx_map(t, u):
    """Returns the pyx snippet for a map of type <t, u>."""
    iterkey = c2py_exprs[t].format(var="deref(inow).first")
    convval = c2py_exprs[u].format(var="v")
    initkey = py2c_exprs[t].format(var="key")
    initval = py2c_exprs[u].format(var="value")
    return _pyxmap.format(tclsname=class_names[t], uclsname=class_names[u],
                          thumname=human_names[t], uhumname=human_names[u],
                          tctype=ctypes[t], uctype=ctypes[u],
                          tpytype=pytypes[t], upytype=pytypes[u],
                          tcytype=cytypes[t], ucytype=cytypes[u],
                          iterkey=iterkey, convval=convval, 
                          initkey=initkey, initval=initval,)

_pxdmap = """# Map{tclsname}{uclsname}
cdef class MapIter{tclsname}{uclsname}(object):
    cdef cpp_map[{tctype}, {uctype}].iterator * iter_now
    cdef cpp_map[{tctype}, {uctype}].iterator * iter_end
    cdef void init(MapIter{tclsname}{uclsname}, cpp_map[{tctype}, {uctype}] *)

cdef class _Map{tclsname}{uclsname}:
    cdef cpp_map[{tctype}, {uctype}] * map_ptr
    cdef public bint _free_map


"""
def genpxd_map(t, u):
    """Returns the pxd snippet for a set of type t."""
    return _pxdmap.format(tclsname=class_names[t], uclsname=class_names[u],
                          thumname=human_names[t], uhumname=human_names[u],
                          tctype=ctypes[t], uctype=ctypes[u],)


_testmap = """# Map{tclsname}{uclsname}
def test_map_{tfncname}_{ufncname}():
    m = conv.Map{tclsname}{uclsname}()
    m[{0}] = {4}
    m[{1}] = {5}
    assert_equal(len(m), 2)
    assert_equal(m[{1}], {5})

    m = conv.Map{tclsname}{uclsname}({{{2}: {6}, {3}: {7}}})
    assert_equal(len(m), 2)
    assert_equal(m[{2}], {6})

    n = conv.Map{tclsname}{uclsname}(m, False)
    assert_equal(len(n), 2)
    assert_equal(n[{2}], {4})

    # points to the same underlying map
    n[{1}] = {5}
    assert_equal(m[{1}], {5})

"""
def gentest_map(t, u):
    """Returns the test snippet for a set of type t."""
    return _testmap.format(*[repr(i) for i in testvals[t] + testvals[u][::-1]], 
                           tclsname=class_names[t], uclsname=class_names[u],
                           tfncname=func_names[t], ufncname=func_names[u])


#
# Python <-> Map Cython Converter Functions
#


_pyxpy2cmap = '''# <{thumname}, {uhumname}> conversions
cdef cpp_map[{tctype}, {uctype}] dict_to_map_{tfncname}_{ufncname}(dict pydict):
    cdef cpp_map[{tctype}, {uctype}] cppmap = cpp_map[{tctype}, {uctype}]()
    for key, value in pydict.items():
        cppmap[{initkey}] = {initval}
    return cppmap

cdef dict map_to_dict_{tfncname}_{ufncname}(cpp_map[{tctype}, {uctype}] cppmap):
    pydict = {{}}
    cdef cpp_map[{tctype}, {uctype}].iterator mapiter = cppmap.begin()
    while mapiter != cppmap.end():
        pydict[{iterkey}] = {iterval}
        inc(mapiter)
    return pydict
'''
def genpyx_py2c_map(t, u):
    """Returns the pyx snippet for a map of type <t, u>."""
    iterkey = c2py_exprs[t].format(var="deref(mapiter).first")
    iterval = c2py_exprs[u].format(var="deref(mapiter).second")
    initkey = py2c_exprs[t].format(var="key")
    initval = py2c_exprs[u].format(var="value")
    return _pyxpy2cmap.format(tclsname=class_names[t], uclsname=class_names[u],
                              thumname=human_names[t], uhumname=human_names[u],
                              tctype=ctypes[t], uctype=ctypes[u],
                              tpytype=pytypes[t], upytype=pytypes[u],
                              tcytype=cytypes[t], ucytype=cytypes[u],
                              iterkey=iterkey, iterval=iterval, 
                              initkey=initkey, initval=initval,
                              tfncname=func_names[t], ufncname=func_names[u],
                              )

_pxdpy2cmap = """# <{thumname}, {uhumname}> conversions
cdef cpp_map[{tctype}, {uctype}] dict_to_map_{tfncname}_{ufncname}(dict)
cdef dict map_to_dict_{tfncname}_{ufncname}(cpp_map[{tctype}, {uctype}])

"""
def genpxd_py2c_map(t, u):
    """Returns the pxd snippet for a set of type t."""
    return _pxdpy2cmap.format(tclsname=class_names[t], uclsname=class_names[u],
                              thumname=human_names[t], uhumname=human_names[u],
                              tctype=ctypes[t], uctype=ctypes[u],
                              tfncname=func_names[t], ufncname=func_names[u])


def gentest_py2c_map(t, u):
    return ""


#
# Controlers 
#

_pyxheader = """###################
###  WARNING!!! ###
###################
# This file has been autogenerated

# Cython imports
from libcpp.utility cimport pair
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.vector cimport vector as cpp_vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libc.stdlib cimport malloc, free

# Python Imports
import collections

cimport numpy as np
import numpy as np

# Local imports
include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
cimport extra_types


cdef extra_types.complex_t _py2c_complex(object pyv):
    cdef extra_types.complex_t cv
    pyv = complex(pyv)
    cv = extra_types.complex_t()
    cv.re = pyv.real
    cv.im = pyv.imag
    return cv


"""
def genpyx(template, header=None):
    """Returns a string of a pyx file representing the given template."""
    pyxfuncs = dict([(k[7:], v) for k, v in globals().items() \
                    if k.startswith('genpyx_') and callable(v)])
    pyx = _pyxheader if header is None else header
    for t in template:
        pyx += pyxfuncs[t[0]](*t[1:]) + "\n\n" 
    return pyx


_pxdheader = """###################
###  WARNING!!! ###
###################
# This file has been autogenerated

# Cython imports
from libcpp.utility cimport pair
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.vector cimport vector as cpp_vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# Python Imports
cimport numpy as np

# Local imports
include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
cimport extra_types


cdef extra_types.complex_t py2c_complex(object)

"""
def genpxd(template, header=None):
    """Returns a string of a pxd file representing the given template."""
    pxdfuncs = dict([(k[7:], v) for k, v in globals().items() \
                    if k.startswith('genpxd_') and callable(v)])
    pxd = _pxdheader if header is None else header
    for t in template:
        pxd += pxdfuncs[t[0]](*t[1:]) + "\n\n" 
    return pxd


_testheader = '''"""Tests the part of stlconverters that is accessible from Python."""
###################
###  WARNING!!! ###
###################
# This file has been autogenerated

from unittest import TestCase
import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

import os
import numpy  as np
import tables as tb

import pyne.stlconverters as conv


'''
def gentest(template, header=None):
    """Returns a string of a test file representing the given template."""
    testfuncs = dict([(k[8:], v) for k, v in globals().items() \
                    if k.startswith('gentest_') and callable(v)])
    test = _testheader if header is None else header
    for t in template:
        test += testfuncs[t[0]](*t[1:]) + "\n\n" 
    return test


def genfiles(template, fname='temp', pxdname=None, testname=None, 
             pyxheader=None, pxdheader=None, testheader=None):
    """Generates all cython source files needed to create the wrapper."""
    # munge some filenames
    fname = fname[:-4] if fname.endswith('.pyx') else fname
    pxdname = fname if pxdname is None else pxdname
    pxdname = pxdname + '.pxd' if not pxdname.endswith('.pxd') else pxdname
    testname = 'test_' + fname if testname is None else testname
    testname = testname + '.py' if not testname.endswith('.py') else testname
    fname += '.pyx'

    pyx = genpyx(template, pyxheader)
    with open(fname, 'w') as f:
        f.write(pyx)

    pxd = genpxd(template, pxdheader)
    with open(pxdname, 'w') as f:
        f.write(pxd)

    test = gentest(template, testheader)
    with open(testname, 'w') as f:
        f.write(test)

if __name__ == "__main__":
    #t = [('set', 'int')]
    #t = [('set', 'str')]
    t = [('py2c_map', 'int', 'int')]
    #print gentest(t)
    #print genpxd(t)
    print genpyx(t)
