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
    t = [('set', 'str')]
    print gentest(t)
