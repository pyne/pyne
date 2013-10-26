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
include "include/cython_version.pxi"
IF CYTHON_VERSION_MAJOR == 0 and CYTHON_VERSION_MINOR >= 17:
    from libcpp.map cimport map as cpp_map
    from libcpp.set cimport set as cpp_set
    from libcpp.string cimport string as std_string
    from libcpp.utility cimport pair as cpp_pair
ELSE:
    from pyne._includes.libcpp.map cimport map as cpp_map
    from pyne._includes.libcpp.set cimport set as cpp_set
    from pyne._includes.libcpp.string cimport string as std_string
    from pyne._includes.libcpp.utility cimport pair as cpp_pair
cimport extra_types

cimport pyne.cpp_pyne
cimport pyne.pyne_config
import pyne.pyne_config

cimport pyne.cpp_nucname
cimport pyne.nucname
import pyne.nucname

cimport cpp_data
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv


#
# atomic_mass functions
#
cdef conv._MapIntDouble atomic_mass_map_proxy = conv.MapIntDouble(False)
atomic_mass_map_proxy.map_ptr = &cpp_data.atomic_mass_map
atomic_mass_map = atomic_mass_map_proxy

def atomic_mass(nuc):
    """Finds the atomic mass of a nuclide in [amu].

    Parameters
    ----------
    nuc : int or str
        Input nuclide.

    Returns
    -------
    mass : float
        Atomic mass of this nuclide [amu].

    Notes
    -----
    If the nuclide is not found, the A-number is returned as a float.
    """
    if isinstance(nuc, int):
        mass = cpp_data.atomic_mass(<int> nuc)
    elif isinstance(nuc, basestring):
        mass = cpp_data.atomic_mass(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return mass



#
# natural_abund functions
#
cdef conv._MapIntDouble natural_abund_map_proxy = conv.MapIntDouble(False)
natural_abund_map_proxy.map_ptr = &cpp_data.natural_abund_map
natural_abund_map = natural_abund_map_proxy

# initialize natural_abund_map
cpp_data.natural_abund(<int>10000000)

abundance_by_z = dict([(i, []) for i in range(1,119)])
for zas, abundance in natural_abund_map.items():
    if 0.0 < abundance < 1.0:
        abundance_by_z[zas/10000000].append((zas, abundance))


def natural_abund(nuc):
    """Finds the natural abundance of a nuclide.

    Parameters
    ----------
    nuc : int or str
        Input nuclide.

    Returns
    -------
    abund : float
        Natural abundance of this nuclide.

    Notes
    -----
    If the nuclide is not found, abundance is 0.
    """
    if isinstance(nuc, int):
        abund = cpp_data.natural_abund(<int> nuc)
    elif isinstance(nuc, basestring):
        abund = cpp_data.natural_abund(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return abund



#
# scattering length functions
#
cdef conv._MapIntComplex b_coherent_map_proxy = conv.MapIntComplex(False)
b_coherent_map_proxy.map_ptr = &cpp_data.b_coherent_map
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



cdef conv._MapIntComplex b_incoherent_map_proxy = conv.MapIntComplex(False)
b_incoherent_map_proxy.map_ptr = &cpp_data.b_incoherent_map
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



cdef conv._MapIntDouble b_map_proxy = conv.MapIntDouble(False)
b_map_proxy.map_ptr = &cpp_data.b_map
b_map = b_map_proxy

def b(nuc):
    """Finds the bound scattering length of a nuclide in [cm].

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    b : float
        Bound scattering length of nuc [cm].

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
cdef conv._MapIntDouble half_life_map_proxy = conv.MapIntDouble(False)
half_life_map_proxy.map_ptr = &cpp_data.half_life_map
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



cdef conv._MapIntDouble decay_const_map_proxy = conv.MapIntDouble(False)
decay_const_map_proxy.map_ptr = &cpp_data.decay_const_map
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



def branch_ratio(from_nuc, to_nuc):
    """Finds a branch ratio for a from -> to nuclide pair [fraction].

    Parameters
    ----------
    from_nuc : int or str 
        Parent nuclide.
    to_nuc : int or str 
        Child nuclide.

    Returns
    -------
    br : float
        Branch ratio of this nuclide pair [fraction].

    Notes
    -----
    If this pair is not found, it is assumed to be impossible, and the branch ratio
    is set to zero.
    """
    if isinstance(from_nuc, int):
        fn = pyne.cpp_nucname.id(<int> from_nuc)
    elif isinstance(from_nuc, basestring):
        fn = pyne.cpp_nucname.id(std_string(<char *> from_nuc))
    else:
        raise pyne.nucname.NucTypeError(from_nuc)

    if isinstance(to_nuc, int):
        tn = pyne.cpp_nucname.id(<int> to_nuc)
    elif isinstance(to_nuc, basestring):
        tn = pyne.cpp_nucname.id(std_string(<char *> to_nuc))
    else:
        raise pyne.nucname.NucTypeError(to_nuc)

    br = cpp_data.branch_ratio(cpp_pair[int, int](fn, tn))
    return br


cdef conv._MapIntDouble state_energy_map_proxy = conv.MapIntDouble(False)
state_energy_map_proxy.map_ptr = &cpp_data.state_energy_map
state_energy_map = state_energy_map_proxy

def state_energy(nuc):
    """Finds the excitation energy [MeV] of a nuclide in a given state.

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    se : float
        Excitation energy of this nuclide [MeV].

    Notes
    -----
    If the nuclide is not found, the nuclide is assumed to be stable.
    """
    if isinstance(nuc, int):
        se = cpp_data.state_energy(<int> nuc)
    elif isinstance(nuc, basestring):
        se = cpp_data.state_energy(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return se


def decay_children(nuc):
    """Finds the decay children of a nuclide.

    Parameters
    ----------
    nuc : int or str 
        Input nuclide.

    Returns
    -------
    dc : set of ints
        Decay children in id form.

    Notes
    -----
    If the nuclide is not found or is stable, the empty set is returned.
    """
    cdef conv._SetInt dc = conv.SetInt()

    if isinstance(nuc, int):
        dc.set_ptr[0] = cpp_data.decay_children(<int> nuc)
    elif isinstance(nuc, basestring):
        dc.set_ptr[0] = cpp_data.decay_children(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return dc
