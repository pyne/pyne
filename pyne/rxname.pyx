"""Python wrapper for rxname library."""

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from libc.string cimport const_char
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.string cimport string as std_string
ELSE:
    from pyne._includes.libcpp.string cimport string as std_string
cimport extra_types
cimport pyne.cpp_pyne
cimport pyne.pyne_config
import pyne.pyne_config

cimport cpp_nucname
cimport cpp_rxname
cimport pyne.stlconverters as conv
import pyne.stlconverters as conv

# names
cdef conv._SetStr names_proxy = conv.SetStr(False)
names_proxy.set_ptr = &cpp_rxname.names
names = names_proxy

# id_name
cdef conv._MapUIntStr id_name_proxy = conv.MapUIntStr(False)
id_name_proxy.map_ptr = &cpp_rxname.id_name
id_name = id_name_proxy

# name_id
cdef conv._MapStrUInt name_id_proxy = conv.MapStrUInt(False)
name_id_proxy.map_ptr = &cpp_rxname.name_id
name_id = name_id_proxy

# id_mt
cdef conv._MapUIntUInt id_mt_proxy = conv.MapUIntUInt(False)
id_mt_proxy.map_ptr = &cpp_rxname.id_mt
id_mt = id_mt_proxy

# mt_id
cdef conv._MapUIntUInt mt_id_proxy = conv.MapUIntUInt(False)
mt_id_proxy.map_ptr = &cpp_rxname.mt_id
mt_id = mt_id_proxy

# labels
cdef conv._MapUIntStr labels_proxy = conv.MapUIntStr(False)
labels_proxy.map_ptr = &cpp_rxname.labels
labels = labels_proxy

# docs
cdef conv._MapUIntStr docs_proxy = conv.MapUIntStr(False)
docs_proxy.map_ptr = &cpp_rxname.docs
docs = docs_proxy



def hash(char * s):
    """Hashes a string to be used as a reaction id.  This hash specifically reserves
    the h < 1000 for MT numbers.
    """
    return int(cpp_rxname.hash(<const_char *> s))


def name(x, y=None, char * z="n"):
    """name(x)

    Gives the unique reaction name.  Note that this name follows the 'natural naming' 
    convention and may be used as a variable in most languages.

    Parameters
    ----------
    x : str, int, or long
        name, abbreviation, id, MT number, or from nuclide.
    y : str, int, or None, optional
        to nuclide.
    z : str, optional
        incident particle type ("n", "p", ...) when x and y are nuclides.

    Returns
    -------
    n : str
        a unique reaction name.
    """
    cdef int from_nuc, to_nuc
    if y is None:
        if isinstance(x, basestring):
            n = cpp_rxname.name(std_string(<char *> x))
        elif isinstance(x, int):
            n = cpp_rxname.name(<extra_types.uint> long(x))
        elif isinstance(x, long):
            n = cpp_rxname.name(<extra_types.uint> x)
    else:
        if isinstance(x, basestring):
            from_nuc = cpp_nucname.zzaaam(std_string(<char *> x))
        elif isinstance(x, int):
            from_nuc = cpp_nucname.zzaaam(<int> x)
        if isinstance(y, basestring):
            to_nuc = cpp_nucname.zzaaam(std_string(<char *> y))
        elif isinstance(y, int):
            to_nuc = cpp_nucname.zzaaam(<int> y)
        n = cpp_rxname.name(from_nuc, to_nuc, std_string(z))
    return n
