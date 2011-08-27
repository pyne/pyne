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
cimport pyne.cpp_pyne
cimport pyne.pyne_config
import pyne.pyne_config

cimport pyne.cpp_nucname
cimport pyne.nucname
import pyne.nucname

cimport cpp_data
cimport pyne.stlconverters as conv
import pyne.stlconverters as conv


#
# nuc_weight Functions
#
cdef conv._MapProxyIntDouble nuc_weight_map_proxy = conv.MapProxyIntDouble()
nuc_weight_map_proxy.init(&cpp_data.nuc_weight_map)
nuc_weight_map = nuc_weight_map_proxy

def nuc_weight(nuc):
    """Finds the weight of a nuclide in [amu].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    weight : float
        Atomic weight of this nuclide [amu].
    """
    if isinstance(nuc, basestring):
        weight = cpp_data.nuc_weight(<char *> nuc)
    elif isinstance(nuc, int):
        weight = cpp_data.nuc_weight(<int> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return weight


