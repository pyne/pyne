"""Python wrapper for nucname library."""
# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libcpp.string cimport string as std_string
#from cython cimport pointer

# Python imports
#from collections import Iterable
from warnings import warn

# local imports
cimport extra_types
cimport pyne.cpp_pyne
cimport pyne.pyne_config
import pyne.pyne_config

cimport cpp_nucname
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv

warn(__name__ + " is not yet V&V compliant.", ImportWarning)
