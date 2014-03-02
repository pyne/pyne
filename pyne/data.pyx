"""Python wrapper for nucname library."""
# Python imports 
#from collections import Iterable

# Cython imports
from libcpp.map cimport map
from libcpp.set cimport set as cpp_set
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string
from libcpp.utility cimport pair as cpp_pair
#from cython cimport pointer

import numpy as np
cimport numpy as np

# local imports 
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

# Mathematical constants
pi = cpp_data.pi
"""Mathematical constant pi."""

N_A = cpp_data.N_A
"""Avogadro constant."""

barns_per_cm2 = cpp_data.barns_per_cm2
"""Barns per centimeter squared."""

cm2_per_barn = cpp_data.cm2_per_barn
"""Centimeter squared per barn."""

sec_per_day = cpp_data.sec_per_day
"""The number of seconds in a canonical day."""


#
# hash map and initialization
#
cdef conv._MapStrStr data_checksums_proxy = conv.MapStrStr(False)
data_checksums_proxy.map_ptr = &cpp_data.data_checksums
data_checksums = data_checksums_proxy


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
cpp_data.natural_abund(<int> 10000000)

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
# Fission product yield data
#


def fpyield(from_nuc, to_nuc, source=0, get_errors=False):
    """Finds the fission product yield for a (parent, child) nuclide pair [fraction].

    Parameters
    ----------
    from_nuc : int or str 
        Parent nuclide.
    to_nuc : int or str 
        Child nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'WIMSD': 0, 'NDS_THERMAL' : 1, 'NDS_FAST' : 2, 'NDS_14MEV' : 3
    get_errors : boolean
        return the error in the value if possible or 0

    Returns
    -------
    fpy : float
        Fractional yield of this nuclide pair [unitless].

    Notes
    -----
    If this pair is not found, it is assumed to be impossible, and the yield
    is set to zero.
    """
    srcmap = {'WIMSD': 0, 'NDS_THERMAL': 1, 'NDS_FAST': 2, 'NDS_14MEV': 3}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 3:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')
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
    fpy = cpp_data.fpyield(cpp_pair[int, int](fn, tn), <int> source, get_errors)
    return fpy


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

def metastable_id(nuc, level=1):
    """
    return the nuc_id of a metastable state

    Parameters
    ----------
    nuc : int
        Input nuclide
    level : int
        integer metastable state
    """
    cpp_data.metastable_id(<int> nuc, <int> level)


gammas_dtype = np.dtype([
    ('energy', float),
    ('energy_err', float),
    ('photon_intensity', float),
    ('photon_intensity_err', float),
    ('conv_intensity', float),
    ('conv_intensity_err', float),
    ('total_intensity', float),
    ('total_intensity_err', float),
    ('from_nuc', int),
    ('to_nuc', int),
    ('parent_nuc', int),
    ('k_conv_e', float),
    ('l_conv_e', float),
    ('m_conv_e', float),
    ])


def get_gammas_by_en(en, pm = 1.0):
    """
    return a list of gamma rays with energies between en + pm and en - pm

    Parameters
    ----------
    en : double
        energy in kev
    pm : double
        neighboring region to search in keV
    """

    cdef cpp_data.gamma_struct * gamma_st
    len = cpp_data.gamma_data_byen(<double> en, <double> pm, gamma_st)
    cdef np.ndarray arr = np.require(np.ndarray(len,dtype=gammas_dtype),requirements=['C','A'])
    for i in range(len):
        arr[i]['energy'] = gamma_st[i].energy
        arr[i]['energy_err'] = gamma_st[i].energy_err
        arr[i]['photon_intensity'] = gamma_st[i].photon_intensity
        arr[i]['photon_intensity_err'] = gamma_st[i].photon_intensity_err
        arr[i]['conv_intensity'] = gamma_st[i].conv_intensity
        arr[i]['conv_intensity_err'] = gamma_st[i].conv_intensity_err
        arr[i]['total_intensity'] = gamma_st[i].total_intensity
        arr[i]['total_intensity_err'] = gamma_st[i].total_intensity_err
        arr[i]['from_nuc'] = gamma_st[i].from_nuc
        arr[i]['to_nuc'] = gamma_st[i].to_nuc
        arr[i]['parent_nuc'] = gamma_st[i].parent_nuc
        arr[i]['k_conv_e'] = gamma_st[i].k_conv_e
        arr[i]['l_conv_e'] = gamma_st[i].l_conv_e
        arr[i]['l_conv_e'] = gamma_st[i].m_conv_e

    return arr

def get_gammas_by_parent(id):
    """
    return a list of gamma rays with a given parent

    Parameters
    ----------
    id : int
        id of parent nuclide
    """

    cdef cpp_data.gamma_struct * gamma_st
    cpp_data.gamma_data_byparent(<int> id, gamma_st)


def get_alphas_by_en(en, pm = 1.0):
    """
    return a list of alpha rays with energies between en + pm and en - pm

    Parameters
    ----------
    en : double
        energy in kev
    pm : double
        neighboring region to search in keV
    """

    cdef cpp_data.alpha_struct * alpha_st
    cpp_data.alpha_data_byen(<double> en, <double> pm, alpha_st)

def get_alphas_by_parent(id):
    """
    return a list of alpha rays with a given parent

    Parameters
    ----------
    id : int
        id of parent nuclide
    """

    cdef cpp_data.alpha_struct * alpha_st
    cpp_data.alpha_data_byparent(<int> id, alpha_st)
