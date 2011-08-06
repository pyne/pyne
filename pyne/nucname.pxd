"""Python wrapper for nucname library."""
# Python imports 
#from collections import Iterable

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
#from cython cimport pointer

# local imports 
cimport std
cimport cpp_nucname
cimport pyne.stlconverters as conv
import pyne.stlconverters as conv

#
# Conversion dictionaries
#

#
# Elemental string sets
#


cdef cpp_set[int] zzaaam_set(object nuc_sequence)

