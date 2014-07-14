"""Python wrapper for nucname library."""

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
#from cython cimport pointer
from libcpp.string cimport string as std_string

# Python imports 

# local imports 
cimport pyne.cpp_utils
cimport pyne.pyne_config
import pyne.pyne_config

cimport cpp_nucname
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv

#
# Conversion dictionaries
#

#
# Elemental string sets
#


cdef cpp_set[int] id_set(object nuc_sequence)
cdef cpp_set[int] zzaaam_set(object nuc_sequence)
cdef cpp_set[int] zzzaaa_set(object nuc_sequence)
