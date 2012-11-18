"""Python wrapper for rxname library."""

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
cimport pyne.cpp_pyne
cimport pyne.pyne_config
import pyne.pyne_config

cimport cpp_rxname
cimport pyne.stlconverters as conv
import pyne.stlconverters as conv

# names
cdef conv._SetStr names_proxy = conv.SetStr(False)
names_proxy.set_ptr = &cpp_rxname.names
names = names_proxy

# labels
cdef conv._MapUIntStr labels_proxy = conv.MapUIntStr(False)
labels_proxy.map_ptr = &cpp_rxname.labels
labels = labels_proxy

# id_name
cdef conv._MapUIntStr id_name_proxy = conv.MapUIntStr(False)
id_name_proxy.map_ptr = &cpp_rxname.id_name
id_name = id_name_proxy

# name_id
cdef conv._MapStrUInt name_id_proxy = conv.MapStrUInt(False)
name_id_proxy.map_ptr = &cpp_rxname.name_id
name_id = name_id_proxy


