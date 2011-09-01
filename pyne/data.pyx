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
cimport extra_types

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
# nuc_weight functions
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

    Notes
    -----
    If the nuclide is not found, the A-number is returned as a float.
    """
    if isinstance(nuc, int):
        weight = cpp_data.nuc_weight(<int> nuc)
    elif isinstance(nuc, basestring):
        weight = cpp_data.nuc_weight(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return weight



#
# scattering length functions
#
cdef conv._MapProxyIntComplex b_coherent_map_proxy = conv.MapProxyIntComplex()
b_coherent_map_proxy.init(&cpp_data.b_coherent_map)
b_coherent_map = b_coherent_map_proxy


def b_coherent(nuc):
    """Finds the coherent bound scattering length of a nuclide in [cm].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    bc : complex
        Coherent bound scattering length of nuc [cm].

    Notes
    -----
    If nuc is not found, the value for a nuclide with the same A-number 
    is used instead. If still no value is found, the an isotope of the 
    same element as nuc is used.  If still no values are found, zero is
    returned.
    """
    cdef extra_types.complex_t value

    if isinstance(nuc, int):
        value = cpp_data.b_coherent(<int> nuc)
    elif isinstance(nuc, basestring):
        value = cpp_data.b_coherent(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return complex(float(value.re), float(value.im))



cdef conv._MapProxyIntComplex b_incoherent_map_proxy = conv.MapProxyIntComplex()
b_incoherent_map_proxy.init(&cpp_data.b_incoherent_map)
b_incoherent_map = b_incoherent_map_proxy

def b_incoherent(nuc):
    """Finds the incoherent bound scattering length of a nuclide in [cm].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    bi : complex
        Incoherent bound scattering length of nuc [cm].

    Notes
    -----
    If nuc is not found, the value for a nuclide with the same A-number 
    is used instead. If still no value is found, the an isotope of the 
    same element as nuc is used.  If still no values are found, zero is
    returned.
    """
    cdef extra_types.complex_t value

    if isinstance(nuc, int):
        value = cpp_data.b_incoherent(<int> nuc)
    elif isinstance(nuc, basestring):
        value = cpp_data.b_incoherent(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return complex(float(value.re), float(value.im))



cdef conv._MapProxyIntDouble b_map_proxy = conv.MapProxyIntDouble()
b_map_proxy.init(&cpp_data.b_map)
b_map = b_map_proxy

def b(nuc):
    """Finds the bound scattering length of a nuclide in [cm].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    bc : float
        Coherent bound scattering length of nuc [cm].

    Notes
    -----
    If nuc is not found, the value for a nuclide with the same A-number 
    is used instead. If still no value is found, the an isotope of the 
    same element as nuc is used.  If still no values are found, zero is
    returned.

    This value is computed from the coherent and incoherent scattering 
    lengths as follows:

    .. math::

        b = \\sqrt{\\left| b_{\mbox{coh}} \\right|^2 + \\left| b_{\mbox{inc}} \\right|^2}

    """
    if isinstance(nuc, int):
        value = cpp_data.b(<int> nuc)
    elif isinstance(nuc, basestring):
        value = cpp_data.b(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return float(value)





#
# decay data functions
#
cdef conv._MapProxyIntDouble half_life_map_proxy = conv.MapProxyIntDouble()
half_life_map_proxy.init(&cpp_data.half_life_map)
half_life_map = half_life_map_proxy

def half_life(nuc):
    """Finds the half-life of a nuclide in [seconds].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    hl : float
        Half-life of this nuclide [seconds].

    Notes
    -----
    If the nuclide is not found, the nuclide is assumed to be stable.
    """
    if isinstance(nuc, int):
        hl = cpp_data.half_life(<int> nuc)
    elif isinstance(nuc, basestring):
        hl = cpp_data.half_life(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return hl



cdef conv._MapProxyIntDouble decay_const_map_proxy = conv.MapProxyIntDouble()
decay_const_map_proxy.init(&cpp_data.decay_const_map)
decay_const_map = decay_const_map_proxy

def decay_const(nuc):
    """Finds the decay constant of a nuclide in [1/seconds].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    dc : float
        Decay constant of this nuclide [1/seconds].

    Notes
    -----
    If the nuclide is not found, the nuclide is assumed to be stable.
    """
    if isinstance(nuc, int):
        dc = cpp_data.decay_const(<int> nuc)
    elif isinstance(nuc, basestring):
        dc = cpp_data.decay_const(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return dc


