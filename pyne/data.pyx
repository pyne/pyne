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
from libc.stdlib cimport free
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
# q_val functions
#
cdef conv._MapIntDouble q_val_map_proxy = conv.MapIntDouble(False)
q_val_map_proxy.map_ptr = &cpp_data.q_val_map
q_val_map = q_val_map_proxy

def q_val(nuc):
    """Finds the Q value of a nuclide in [MeV/fission].

    Parameters
    ----------
    nuc : int or str
        Input nuclide.

    Returns
    -------
    q_val : double
        Q value of this nuclide [MeV/fission].

    Notes
    -----
    If the nuclide is not found, 0 is returned.
    """
    if isinstance(nuc, int):
        q_val = cpp_data.q_val(<int> nuc)
    elif isinstance(nuc, basestring):
        q_val = cpp_data.q_val(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return q_val

# gamma_frac
cdef conv._MapIntDouble gamma_frac_map_proxy = conv.MapIntDouble(False)
gamma_frac_map_proxy.map_ptr = &cpp_data.gamma_frac_map
gamma_frac_map = gamma_frac_map_proxy

def gamma_frac(nuc):
    """Finds the fraction of Q that comes from gammas of a nuclide.

    Parameters
    ----------
    nuc : int or str
        Input nuclide.

    Returns
    -------
    gamma_frac : double
        Fraction of Q that comes from gammas of this nuclide.

    Notes
    -----
    If the nuclide is not found, gamma_frac is 0.
    """
    if isinstance(nuc, int):
        gamma_frac = cpp_data.gamma_frac(<int> nuc)
    elif isinstance(nuc, basestring):
        gamma_frac = cpp_data.gamma_frac(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return gamma_frac


#
# Decay Factor data
#

# external air dose
cdef conv._MapIntDouble ext_air_dose_map_proxy = conv.MapIntDouble(False)
ext_air_dose_map_proxy.map_ptr = &cpp_data.ext_air_dose_map
ext_air_dose_map = ext_air_dose_map_proxy

def ext_air_dose(nuc, source=0):
    """Finds the external air dose factor for a tracked nuclide.

    Parameters
    ----------
    nuc : int or str 
        Parent nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'EPA': 0, 'DOE' : 1, 'GENII' : 2

    Returns
    -------
    ext_air_dose : float
        Dose factor from external air exposure [mrem/hr per Ci/m^3]

    Notes
    -----
    If the nuclide is not found, and the dose factor is set to zero.
    """
    srcmap = {'EPA': 0, 'DOE': 1, 'GENII': 2}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 2:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')

    if isinstance(nuc, int):
        ext_air_dose = cpp_data.ext_air_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        ext_air_dose = cpp_data.ext_air_dose(<char *> nuc, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return ext_air_dose

# ratio
cdef conv._MapIntDouble ratio_map_proxy = conv.MapIntDouble(False)
ratio_map_proxy.map_ptr = &cpp_data.ratio_map
ratio_map = ratio_map_proxy

def ratio(nuc, source=0):
    """Finds ratio of dose from external air to dose from inhalation for a tracked nuclide.

    Parameters
    ----------
    nuc : int or str 
        Parent nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'EPA': 0, 'DOE' : 1, 'GENII' : 2

    Returns
    -------
    ratio : float
        Fraction of dose from external air to dose from inhalation.

    Notes
    -----
    If the nuclide is not found, and the ratio is set to zero.
    """
    srcmap = {'EPA': 0, 'DOE': 1, 'GENII': 2}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 2:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')

    if isinstance(nuc, int):
        ratio = cpp_data.ratio(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        ratio = cpp_data.ratio(<char *> nuc, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return ratio

# external soil dose
cdef conv._MapIntDouble ext_soil_dose_map_proxy = conv.MapIntDouble(False)
ext_soil_dose_map_proxy.map_ptr = &cpp_data.ext_soil_dose_map
ext_soil_dose_map = ext_soil_dose_map_proxy

def ext_soil_dose(nuc, source=0):
    """Finds the external soil dose factor for a tracked nuclide.

    Parameters
    ----------
    nuc : int or str 
        Parent nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'EPA': 0, 'DOE' : 1, 'GENII' : 2

    Returns
    -------
    ext_soil_dose : float
        Dose factor from 15 cm of external soil exposure [mrem/hr per Ci/m^2]

    Notes
    -----
    If the nuclide is not found, and the dose factor is set to zero.
    """
    srcmap = {'EPA': 0, 'DOE': 1, 'GENII': 2}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 2:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')

    if isinstance(nuc, int):
        ext_soil_dose = cpp_data.ext_soil_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        ext_soil_dose = cpp_data.ext_soil_dose(<char *> nuc, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return ext_soil_dose
    
# ingestion dose
cdef conv._MapIntDouble ingest_dose_map_proxy = conv.MapIntDouble(False)
ingest_dose_map_proxy.map_ptr = &cpp_data.ingest_dose_map
ingest_dose_map = ingest_dose_map_proxy

def ingest_dose(nuc, source=0):
    """Finds the dose factor due to ingestion for a tracked nuclide.

    Parameters
    ----------
    nuc : int or str 
        Parent nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'EPA': 0, 'DOE' : 1, 'GENII' : 2

    Returns
    -------
    ingest : float
        Dose factor from exposure due to ingestion [mrem/pCi]

    Notes
    -----
    If the nuclide is not found, and the dose factor is set to zero.
    """
    srcmap = {'EPA': 0, 'DOE': 1, 'GENII': 2}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 2:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')

    if isinstance(nuc, int):
        ingest_dose = cpp_data.ingest_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        ingest_dose = cpp_data.ingest_dose(<char *> nuc, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return ingest_dose
    
# fluid_frac
cdef conv._MapIntDouble fluid_frac_map_proxy = conv.MapIntDouble(False)
fluid_frac_map_proxy.map_ptr = &cpp_data.fluid_frac_map
fluid_frac_map = fluid_frac_map_proxy

def fluid_frac(nuc, source=0):
    """Finds fraction of activity that is absorbed by body fluids for a tracked nuclide.

    Parameters
    ----------
    nuc : int or str 
        Parent nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'EPA': 0, 'DOE' : 1, 'GENII' : 2

    Returns
    -------
    fluid_frac : float
        Fraction of activity that is absorbed by body fluids.

    Notes
    -----
    If the nuclide is not found, and the fluid fraction is set to zero.
    """
    srcmap = {'EPA': 0, 'DOE': 1, 'GENII': 2}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 2:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')

    if isinstance(nuc, int):
        fluid_frac = cpp_data.fluid_frac(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        fluid_frac = cpp_data.fluid_frac(<char *> nuc, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return fluid_frac

# inhalation dose
cdef conv._MapIntDouble inhale_dose_map_proxy = conv.MapIntDouble(False)
inhale_dose_map_proxy.map_ptr = &cpp_data.inhale_dose_map
inhale_dose_map = inhale_dose_map_proxy

def inhale_dose(nuc, source=0):
    """Finds the dose factor due to inhalation for a tracked nuclide.

    Parameters
    ----------
    nuc : int or str 
        Parent nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'EPA': 0, 'DOE' : 1, 'GENII' : 2

    Returns
    -------
    inhale_dose : float
        Dose factor from exposure due to inhalation [mrem/pCi]

    Notes
    -----
    If the nuclide is not found, and the dose factor is set to zero.
    """
    srcmap = {'EPA': 0, 'DOE': 1, 'GENII': 2}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 2:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')

    if isinstance(nuc, int):
        inhale_dose = cpp_data.inhale_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        inhale_dose = cpp_data.inhale_dose(<char *> nuc, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return inhale_dose

# lung model
cdef conv._MapIntStr lung_mod_map_proxy = conv.MapIntStr(False)
lung_mod_map_proxy.map_ptr = &cpp_data.lung_mod_map
lung_mod_map = lung_mod_map_proxy

def lung_mod(nuc, source=0):
    """Finds the lung model for the inhalation dose factor for a tracked nuclide.

    Parameters
    ----------
    nuc : int or str 
        Parent nuclide.
    source : int or str
        The int or corresponding dictionary key for the source dataset.
        Allowed values are:
        'EPA': 0, 'DOE' : 1, 'GENII' : 2

    Returns
    -------
    lung_mod : string
        Model of lung used for calculation (D (days), W (weeks), or Y (years)).

    Notes
    -----
    If the nuclide is not found, and the lung model is set to zero.
    """
    srcmap = {'EPA': 0, 'DOE': 1, 'GENII': 2}
    if isinstance(source, str):
        sourceint = srcmap[source]
    elif isinstance(source, int):
        if 0 <= source <= 2:
            sourceint = source
        else:
            raise ValueError
    else:
        raise ValueError('Only ints or strings are accepted')

    if isinstance(nuc, int):
        lung_mod = cpp_data.lung_mod(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        lung_mod = cpp_data.lung_mod(<char *> nuc, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return lung_mod


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

    Returns
    -------
    nuc : int
        nuc_id of metastable state
    """
    return cpp_data.metastable_id(<int> nuc, <int> level)

def decay_half_life(from_nuc, to_nuc):
    """
    Returns the half life from ENSDF decay dataset parent data

    Parameters
    ----------
    from_nuc : int
        parent nuclide
    to_nuc : int
        child nuclide

    Returns
    -------
    half_life : double
        half life in seconds
    error : double
        Error in seconds
    """
    half_life, error = cpp_data.decay_half_life(cpp_pair[int,int](from_nuc, to_nuc)) 
    return half_life, error

def decay_half_life_byparent(parent):
    """
    Returns a list half lives from ENSDF decay dataset parent data

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    half_lives : array of pairs
        An array of half lives and their errors for a given parent nuclide
    """
    half_lives = cpp_data.decay_half_lifes(<int> parent)
    return half_lives

def decay_branch_ratio(from_nuc, to_nuc):
    """
    Returns the branch ratio from ENSDF decay dataset data

    Parameters
    ----------
    from_nuc : int
        parent nuclide
    to_nuc : int
        child nuclide

    Returns
    -------
    ratio : double
        branching ratio
    """
    ratio = cpp_data.decay_branch_ratio(cpp_pair[int,int](from_nuc, to_nuc)) 
    return ratio

def decay_branch_ratio_byparent(parent):
    """
    Returns a list branching ratios from ENSDF decay dataset data for a given
    parent.

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of doubles
        An array of half lives and their errors for a given parent nuclide
    """
    ratios = cpp_data.decay_branch_ratios(<int> parent)
    return ratios

def decay_photon_branch_ratio(from_nuc, to_nuc):
    """
    Returns the photon branch ratio from ENSDF decay dataset data

    Parameters
    ----------
    from_nuc : int
        parent nuclide
    to_nuc : int
        child nuclide

    Returns
    -------
    ratio : double
        photon branching ratio
    """
    ratio, error = \
    cpp_data.decay_photon_branch_ratio(cpp_pair[int,int](from_nuc, to_nuc)) 
    return ratio, error

def decay_photon_branch_ratio_byparent(parent):
    """
    Returns a list of photon branch ratios from ENSDF decay dataset parent data

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    half_lives : array of pairs
        An array of photon branching ratios and their errors for a given parent
        nuclide
    """
    arr = cpp_data.decay_photon_branch_ratios(<int> parent)
    return arr

def decay_beta_branch_ratio(from_nuc, to_nuc):
    """
    Returns the  branch ratio from ENSDF decay dataset data

    Parameters
    ----------
    from_nuc : int
        parent nuclide
    to_nuc : int
        child nuclide

    Returns
    -------
    ratio : double
         branching ratio
    """
    ratio, error = \
    cpp_data.decay_beta_branch_ratio(cpp_pair[int,int](from_nuc, to_nuc)) 
    return ratio, error

def decay_beta_branch_ratio_byparent(parent):
    """
    Returns a list of beta branch ratios from ENSDF decay dataset parent data

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of beta branching ratios and their errors for a given parent
        nuclide
    """
    ratios = cpp_data.decay_beta_branch_ratios(<int> parent)
    return ratios

def gamma_energy(parent):
    """
    Returns a list of gamma ray energies from ENSDF decay dataset from a given 
    parent

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray energies and errors
    """
    return cpp_data.gamma_energy(<int> parent)

def gamma_photon_intensity(parent):
    """
    Returns a list of gamma ray photon intensities from ENSDF decay dataset 
    from a given parent

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray photon intensities and errors
    """
    return cpp_data.gamma_photon_intensity(<int> parent)
    
def gamma_conversion_intensity(parent):
    """
    Returns a list of gamma ray conversion intensities from ENSDF decay dataset
    from a given parent

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray conversion intensities and errors
    """
    return cpp_data.gamma_conversion_intensity(<int> parent)
    
def gamma_total_intensity(parent):
    """
    Returns a list of gamma ray total intensities from ENSDF decay dataset from
    a given parent

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray total intensities and errors
    """
    return cpp_data.gamma_total_intensity(<int> parent)

def gamma_from_to_byparent(parent):
    """
    Returns a list of gamma ray level pairs from ENSDF decay dataset from a
    given parent. This makes it possible to calculate coincidence rates.

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray level pairs in nuc_id form
    """
    return cpp_data.gamma_from_to(<int> parent)
    
def gamma_from_to_byen(en, enerror=None):
    """
    Returns a list of gamma ray level pairs from ENSDF decay dataset 
    based on gamma-ray energy. 

    Parameters
    ----------
    en : double
        gamma ray energy in keV
    enerror : double
        gamma ray energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray level pairs in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.gamma_from_to(<double> en,<double> enerror)   

def gamma_parent(en, enerror=None):
    """
    Returns a list of gamma ray parents from ENSDF decay dataset 
    based on gamma-ray energy. 

    Parameters
    ----------
    en : double
        gamma ray energy in keV
    enerror : double
        gamma ray energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of ints
        An array of gamma ray parents in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.gamma_parent(<double> en, <double> enerror)


def alpha_energy(parent):
    """
    Returns a list of alpha energies from ENSDF decay dataset from a given 
    parent

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of alpha energies and errors
    """
    return cpp_data.alpha_energy(<int> parent)

def alpha_intensity(parent):
    """
    Returns a list of alpha intensities from ENSDF decay dataset from a given 
    parent

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of alpha intensities and errors
    """
    return cpp_data.alpha_intensity(<int> parent) 

def alpha_parent(en, enerror=None):
    """
    Returns a list of alpha parents from ENSDF decay dataset 
    based on alpha energy. 

    Parameters
    ----------
    en : double
        alpha energy in keV
    enerror : double
        alpha energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of ints
        An array of alpha parents in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.alpha_parent(<double> en, <double> enerror)

def alpha_child_byen(en, enerror=None):
    """
    Returns a list of alpha children from ENSDF decay dataset 
    based on alpha energy.

    Parameters
    ----------
    en : double
        alpha energy in keV
    enerror : double
        alpha energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of ints
        An array of alpha children in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.alpha_child(<double> en, <double> enerror)
    
def alpha_child_byparent(parent):
    """
    Returns a list of alpha children from ENSDF decay dataset 
    based on alpha parent.

    Parameters
    ----------
    parent : int
        parent nuclide in nuc_id form

    Returns
    -------
    ratios : array of ints
        An array of alpha children in nuc_id form
    """
    return cpp_data.alpha_child(<int> parent)

def beta_endpoint_energy(parent):
    """
    Returns a list of beta endpoint energies from ENSDF decay dataset 
    based on parent nuclide.

    Parameters
    ----------
    parent : int
        parent nuclide in nuc_id form

    Returns
    -------
    ratios : array of ints
        An array of beta endpoint energies and errors
    """
    return cpp_data.beta_endpoint_energy(<int> parent)

def beta_average_energy(parent):
    """
    Returns a list of beta average energies from ENSDF decay dataset 
    based on parent nuclide.

    Parameters
    ----------
    parent : int
        parent nuclide in nuc_id form

    Returns
    -------
    ratios : array of ints
        An array of beta average energies and errors
    """
    return cpp_data.beta_average_energy(<int> parent)

def beta_intensity(parent):
    """
    Returns a list of beta intensities from ENSDF decay dataset 
    based on parent nuclide.

    Parameters
    ----------
    parent : int
        parent nuclide in nuc_id form

    Returns
    -------
    ratios : array of ints
        An array of beta intensities and errors
    """
    return cpp_data.beta_intensity(<int> parent) 

def beta_parent(en, enerror=None):
    """
    Returns a list of beta minus parents from ENSDF decay dataset 
    based on beta energy.

    Parameters
    ----------
    en : double
        beta- energy in keV
    enerror : double
        beta- energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of ints
        An array of beta minus parents in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.beta_parent(<double> en, <double> enerror)

def beta_child_byen(en, enerror=None):
    """
    Returns a list of beta minus children from ENSDF decay dataset 
    based on beta energy.

    Parameters
    ----------
    en : double
        beta- energy in keV
    enerror : double
        beta- energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of ints
        An array of beta minus children in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.beta_child(<double> en, <double> enerror)
    
def beta_child_byparent(parent):
    """
    Returns a list of beta minus children from ENSDF decay dataset 
    based on parent.

    Parameters
    ----------
    parent : int
        parent nuclide in nuc_id form

    Returns
    -------
    ratios : array of ints
        An array of beta- children in nuc_id form
    """
    return cpp_data.beta_child(<int> parent)

def ecbp_endpoint_energy(parent):
    """
    Returns a list of beta plus endpoint energies from ENSDF decay dataset from
    a given parent.

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of beta plus endpoint energies and errors
    """
    return cpp_data.ecbp_endpoint_energy(<int> parent)

def ecbp_average_energy(parent):
    """
    Returns a list of beta plus average energies from ENSDF decay dataset from
    a given parent.

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of beta plus average energies and errors
    """
    return cpp_data.ecbp_average_energy(<int> parent)

def ec_intensity(parent):
    """
    Returns a list of electron capture intensisities from ENSDF decay dataset
    from a given parent.

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of electron capture intensisities and errors
    """
    return cpp_data.ec_intensity(<int> parent)

def beta_plus_intensity(parent):
    """
    Returns a list of beta plus intensities from ENSDF decay dataset from
    a given parent.

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of pairs
        An array of beta plus intensities and errors
    """
    return cpp_data.bp_intensity(<int> parent) 

def ecbp_parent(en, enerror=None):
    """
    Returns a list of beta plus/electron capture parents from ENSDF decay 
    dataset based on beta energy.

    Parameters
    ----------
    en : double
        beta- energy in keV
    enerror : double
        beta- energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of ints
        An array of beta plus/electron capture children in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.ecbp_parent(<double> en, <double> enerror)

def ecbp_child_byen(en, enerror=None):
    """
    Returns a list of beta plus/electron capture parents from ENSDF decay 
    dataset based on beta energy.

    Parameters
    ----------
    en : double
        beta- energy in keV
    enerror : double
        beta- energy error (range which you want to search) this defaults
        to 1% of the energy if it is not provided

    Returns
    -------
    ratios : array of ints
        An array of beta plus/electron capture children in nuc_id form
    """
    if enerror == None:
        enerror = en * 0.01
    return cpp_data.ecbp_child(<double> en, <double> enerror)
    
def ecbp_child_byparent(parent):
    """
    Returns a list of beta plus children from ENSDF decay dataset 
    based on parent.

    Parameters
    ----------
    parent : int
        parent nuclide in nuc_id form

    Returns
    -------
    ratios : array of ints
        An array of beta+ children in nuc_id form
    """
    return cpp_data.ecbp_child(<int> parent)
