"""Python wrapper for nucname library."""
from __future__ import division, unicode_literals

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

#Standard lib import
from pyne.utils import QA_warn

cimport numpy as np
import numpy as np

# local imports
cimport extra_types

cimport pyne.cpp_utils
cimport pyne.pyne_config
import pyne.pyne_config

cimport pyne.cpp_nucname
cimport pyne.nucname
import pyne.nucname

cimport cpp_data
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv

QA_warn(__name__)

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

MeV_per_K = cpp_data.MeV_per_K
"""Megaelectronvolts per Kelvin."""

MeV_per_MJ = cpp_data.MeV_per_MJ
"""Megaelectronvolts per megajoule."""

Bq_per_Ci = cpp_data.Bq_per_Ci
"""Becquerel per Curie"""

Ci_per_Bq = cpp_data.Ci_per_Bq
"""Curies per Becquerel"""

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
        nuc_bytes = nuc.encode()
        mass = cpp_data.atomic_mass(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
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

abundance_by_z = dict([(i, []) for i in range(1,119)])
for zas, abundance in natural_abund_map.items():
    if 0.0 < abundance < 1.0:
        abundance_by_z[zas//10000000].append((zas, abundance))


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
        abund = cpp_data.natural_abund(<int> pyne.nucname.id(nuc))
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        abund = cpp_data.natural_abund(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
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
        nuc_bytes = nuc.encode()
        q_val = cpp_data.q_val(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
        q_val = cpp_data.q_val(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)
    return q_val


#
# gamma_frac functions
#
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
        nuc_bytes = nuc.encode()
        gamma_frac = cpp_data.gamma_frac(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
        gamma_frac = cpp_data.gamma_frac(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return gamma_frac


#
# simple_xs functions
#

def simple_xs(nuc, rx, energy):
    """Finds the cross section for the given nuclide and reaction in [barns].
    Uses the simple_xs dataset.

    Parameters
    ----------
    nuc : int or str
        Input nuclide.

    rx : int or str
        Input reaction.

    energy : str
        Energy group for reaction.  Must be one of: "thermal",
        "thermal_maxwell_ave", "resonance_integral", "fourteen_MeV",
        "fission_spectrum_ave".

    Returns
    -------
    xs : double
        cross section value for this nuclide and reaction [barns].

    Notes
    -----
    If the nuclide is not found, 0 is returned.
    """
    #make sure we start with unicode strings
    if isinstance(nuc, bytes):
        nuc = nuc.decode("utf-8")
    if isinstance(rx, bytes):
        rx = rx.decode("utf-8")

    if not isinstance(energy, basestring):
        raise ValueError('energy must be string')
    elif not isinstance(nuc, int) and not isinstance(nuc, basestring):
        raise ValueError('nuc must be int or string')
    elif not isinstance(rx, int) and not isinstance(rx, basestring):
        raise ValueError('rx must be int or string')

    energy_bytes = energy.encode()
    if isinstance(nuc, int) and isinstance(rx, int):
        xs = cpp_data.simple_xs(<int> nuc, <int> rx,
                                std_string(<char *> energy_bytes))
    elif isinstance(nuc, int) and isinstance(rx, basestring):
        rxin_bytes = rx.encode()
        xs = cpp_data.simple_xs(<int> nuc, std_string(<char *> rxin_bytes),
                                std_string(<char *> energy_bytes))
    elif isinstance(nuc, basestring) and isinstance(rx, int):
        nucin_bytes = nuc.encode()
        xs = cpp_data.simple_xs(std_string(<char *> nucin_bytes),
                                <int> rx, std_string(<char *> energy_bytes))
    elif isinstance(nuc, basestring) and isinstance(rx, basestring):
        rxin_bytes = rx.encode()
        nucin_bytes = nuc.encode()
        xs = cpp_data.simple_xs(std_string(<char *> nucin_bytes),
                                std_string(<char *> rxin_bytes),
                                std_string(<char *> energy_bytes))

    return xs


#
# Decay Factor data
#

# external air dose

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
    ----
    The only source that provides this data is EPA; all other
    sources will give a value of -1.
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
        raise ValueError('Only ints are accepted')

    if isinstance(nuc, int):
        ext_air_dose = cpp_data.ext_air_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        ext_air_dose = cpp_data.ext_air_dose(<char *> nuc_bytes, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    if ext_air_dose < 0:
        return float('nan')

    return ext_air_dose

# ratio

def dose_ratio(nuc, source=0):
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
    The only source that provides this data is EPA; all other
    sources will give a value of -1.
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
        raise ValueError('Only ints are accepted')

    if isinstance(nuc, int):
        ratio = cpp_data.dose_ratio(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        ratio = cpp_data.dose_ratio(<char *> nuc_bytes, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    if ratio < 0:
        return float('nan')

    return ratio

# external soil dose

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
    If the nuclide is not found, a value of -1 is returned.
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
        raise ValueError('Only ints are accepted')

    if isinstance(nuc, int):
        ext_soil_dose = cpp_data.ext_soil_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        ext_soil_dose = cpp_data.ext_soil_dose(<char *> nuc_bytes, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    if ext_soil_dose < 0:
        return float('nan')

    return ext_soil_dose

# ingestion dose

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
    If the nuclide is not found, a value of -1 is returned.
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
        raise ValueError('Only ints are accepted')

    if isinstance(nuc, int):
        ingest_dose = cpp_data.ingest_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        ingest_dose = cpp_data.ingest_dose(<char *> nuc_bytes, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    if ingest_dose < 0:
        return float('nan')

    return ingest_dose

# fluid_frac

def dose_fluid_frac(nuc, source=0):
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
    If the nuclide is not found, a value of -1 is returned.
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
        raise ValueError('Only ints are accepted')

    if isinstance(nuc, int):
        fluid_frac = cpp_data.dose_fluid_frac(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        fluid_frac = cpp_data.dose_fluid_frac(<char *> nuc_bytes, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    if dose_fluid_frac < 0:
        return float('nan')

    return fluid_frac

# inhalation dose

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
    If the nuclide is not found, a value of -1 is returned.
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
        raise ValueError('Only ints are accepted')

    if isinstance(nuc, int):
        inhale_dose = cpp_data.inhale_dose(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        inhale_dose = cpp_data.inhale_dose(<char *> nuc_bytes, <int> source)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    if inhale_dose < 0:
        return float('nan')

    return inhale_dose

# lung model

def dose_lung_model(nuc, source=0):
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
    If the nuclide is not found, a string of 'Nada' is returned.
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
        raise ValueError('Only ints are accepted')

    if isinstance(nuc, int):
        lung_mod = cpp_data.dose_lung_model(<int> nuc, <int> source)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        lung_mod = cpp_data.dose_lung_model(<char *> nuc_bytes, <int> source)
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
        nuc_bytes = nuc.encode()
        value = cpp_data.b_coherent(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
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
        nuc_bytes = nuc.encode()
        value = cpp_data.b_incoherent(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
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
        nuc_bytes = nuc.encode()
        value = cpp_data.b(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
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
        from_nuc_bytes = from_nuc.encode()
        fn = pyne.cpp_nucname.id(std_string(<char *> from_nuc_bytes))
    elif isinstance(from_nuc, bytes):
        fn = pyne.cpp_nucname.id(std_string(<char *> from_nuc))
    else:
        raise pyne.nucname.NucTypeError(from_nuc)

    if isinstance(to_nuc, int):
        tn = pyne.cpp_nucname.id(<int> to_nuc)
    elif isinstance(to_nuc, basestring):
        to_nuc_bytes = to_nuc.encode()
        tn = pyne.cpp_nucname.id(std_string(<char *> to_nuc_bytes))
    elif isinstance(to_nuc, bytes):
        tn = pyne.cpp_nucname.id(std_string(<char *> to_nuc))
    else:
        raise pyne.nucname.NucTypeError(to_nuc)
    fpy = cpp_data.fpyield(cpp_pair[int, int](fn, tn), <int> source, get_errors)
    return fpy

#
# atomic data functions
#

def calculate_xray_data(nuc, k_conv, l_conv):
    """Calculates X-ray intensities for a given atom with
    k and l conversion intensities

    Parameters
    ----------
    nuc : int or str
        Input nuclide.
    k_conv : float
        k electron converion coefficient arbitrary units
    l_conv : float
        l electron converion coefficient arbitrary units

    Returns
    -------
    arr : vector of pairs
        Vector of pairs containing the four primary X-rays and their
        intensities: Ka1, Ka2, Kb, L
    """
    z = pyne.nucname.znum(nuc)
    return cpp_data.calculate_xray_data(<int> z, <double> k_conv,
                                        <double> l_conv)
#
# decay data functions
#


def half_life(nuc, use_metastable=True):
    """Finds the half-life of a nuclide in [seconds].

    Parameters
    ----------
    nuc : int or str
        Input nuclide, if metastable is false this uses state_id
    use_metastable : bool
        Assume state of input nuc_id refers to metastable state. Defaults to
        True.

    Returns
    -------
    hl : float
        Half-life of this nuclide [seconds].

    Notes
    -----
    If the nuclide is not found, returns nan.
    """
    if use_metastable is True:
        nuc = pyne.nucname.id(nuc)
        ms = nuc % 10000
        nuc = metastable_id(nuc, ms)
        if nuc is None:
            return float('nan')  # nuclide doesn't exist
    if isinstance(nuc, int):
        hl = cpp_data.half_life(<int> nuc)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        hl = cpp_data.half_life(<char *> nuc_bytes)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return hl




def decay_const(nuc, use_metastable=True):
    """Finds the decay constant of a nuclide in [1/seconds].

    Parameters
    ----------
    nuc : int or str
        Input nuclide, if metastable is false this uses state_id
    use_metastable : bool
        Assume state of input nuc_id refers to metastable state. Defaults to
        True.

    Returns
    -------
    dc : float
        Decay constant of this nuclide [1/seconds].

    Notes
    -----
    If the nuclide is not found, the nuclide is assumed to be stable.
    """
    if use_metastable is True:
        nuc = pyne.nucname.id(nuc)
        ms = nuc % 10000
        nuc = metastable_id(nuc, ms)
        if nuc is None:
            return 0.0  # nuclide doesn't exist, assume stable
    if isinstance(nuc, int):
        dc = cpp_data.decay_const(<int> nuc)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        dc = cpp_data.decay_const(<char *> nuc_bytes)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return dc


def branch_ratio(from_nuc, to_nuc, use_metastable=True):
    """Finds a branch ratio for a from -> to nuclide pair [fraction].

    Parameters
    ----------
    from_nuc : int or str
        Parent nuclide, if metastable is false this uses state id
    to_nuc : int or str
        Child nuclide, if metastable is false this uses state id
    use_metastable : bool
        Assume state of input nuc_id refers to metastable state. Defaults to
        True.

    Returns
    -------
    br : float
        Branch ratio of this nuclide pair [fraction].

    Notes
    -----
    If this pair is not found, it is assumed to be impossible, and the branch ratio
    is set to zero.
    """
    if use_metastable is True:
        from_nuc = pyne.nucname.id(from_nuc)
        to_nuc = pyne.nucname.id(to_nuc)
        ms = from_nuc % 10000
        from_nuc = metastable_id(from_nuc, ms)
        if from_nuc is None:
            return 0.0
        ms = to_nuc % 10000
        to_nuc = metastable_id(to_nuc, ms)
        if to_nuc is None:
            return 0.0
    if isinstance(from_nuc, int):
        fn = pyne.cpp_nucname.id(<int> from_nuc)
    elif isinstance(from_nuc, basestring):
        from_nuc_bytes = from_nuc.encode()
        fn = pyne.cpp_nucname.id(std_string(<char *> from_nuc_bytes))
    else:
        raise pyne.nucname.NucTypeError(from_nuc)

    if isinstance(to_nuc, int):
        tn = pyne.cpp_nucname.id(<int> to_nuc)
    elif isinstance(to_nuc, basestring):
        to_nuc_bytes = to_nuc.encode()
        tn = pyne.cpp_nucname.id(std_string(<char *> to_nuc_bytes))
    else:
        raise pyne.nucname.NucTypeError(to_nuc)

    br = cpp_data.branch_ratio(cpp_pair[int, int](fn, tn))
    return br


def state_energy(nuc, use_metastable=True):
    """Finds the excitation energy [MeV] of a nuclide in a given state.

    Parameters
    ----------
    nuc : int or str
        Input nuclide, if metastable is false this uses state id
    use_metastable : bool
        Assume state of input nuc_id refers to metastable state. Defaults to
        True

    Returns
    -------
    se : float
        Excitation energy of this nuclide [MeV].

    Notes
    -----
    If the nuclide is not found, the nuclide is assumed to be stable.
    """
    if use_metastable is True:
        nuc = pyne.nucname.id(nuc)
        ms = nuc % 10000
        nuc = metastable_id(nuc, ms)
        if nuc is None:
            return 0.0  # nuclide doesn't exist, assume stable
    if isinstance(nuc, int):
        se = cpp_data.state_energy(<int> nuc)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        se = cpp_data.state_energy(<char *> nuc_bytes)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return se


def decay_children(nuc, use_metastable=True):
    """Finds the decay children of a nuclide.

    Parameters
    ----------
    nuc : int or str
        Input nuclide, if metastable is false this uses state id
    use_metastable : bool
        Assume state of input nuc_id refers to metastable state. Defaults to
        True

    Returns
    -------
    dc : set of ints
        Decay children in id form.

    Notes
    -----
    If the nuclide is not found or is stable, the empty set is returned.
    """
    if use_metastable is True:
        nuc = pyne.nucname.id(nuc)
        ms = nuc % 10000
        nuc = metastable_id(nuc, ms)
        if nuc is None:
            return set()
    cdef conv._SetInt dc = conv.SetInt()

    if isinstance(nuc, int):
        dc.set_ptr[0] = cpp_data.decay_children(<int> nuc)
    elif isinstance(nuc, basestring):
        nuc_bytes = nuc.encode()
        dc.set_ptr[0] = cpp_data.decay_children(<char *> nuc_bytes)
    elif isinstance(nuc, bytes):
        dc.set_ptr[0] = cpp_data.decay_children(<char *> nuc)
    else:
        raise pyne.nucname.NucTypeError(nuc)

    return dc

def all_children(nuc):
    """
    returns child nuclides from both level and decay data

    Parameters
    ----------
    nuc : int
        input nuclide in state id form

    Returns
    -------
    ids : list of ints
        list of all known decay children from ensdf data.
    """

    return set(decay_children(nuc, False)) | set(decay_data_children(nuc))

def all_branch_ratio(from_nuc, to_nuc):
    """
    returns the branching ratio if its in either level or decay data

    Parameters
    ----------
    from_nuc : int or str
        Parent nuclide, this uses state id
    to_nuc : int or str
        Child nuclide, this uses state id

    Returns
    -------
    br : float
        Branch ratio of this nuclide pair [fraction].
    """
    br1 = branch_ratio(from_nuc, to_nuc, use_metastable=False)
    br2 = decay_branch_ratio(from_nuc, to_nuc)[0]
    if np.isnan(br1) and not np.isnan(br2):
        return br2
    elif not np.isnan(br1) and np.isnan(br2):
        return br1
    elif br2 == 0 and br1 != 0:
        return br1
    elif br1 == 0 and br2 != 0:
        return br2
    else:
        return br2

def id_from_level(nuc, level, special=""):
    """
    return the state_id for input energy level

    Parameters
    ----------
    nuc : int
        Input nuclide
    level : double
        energy level of state
    special : str
        special level denotation. This is a single A-Z character corresponding
        to a group of levels and associated gammas with no reference to the GS
    Returns
    -------
    nuc : int
        state_id of state
    """
    cdef std_string spc
    if len(special) == 1:
        spc = special[0].encode('UTF-8')
    if level is not None and (level > 0.0 or len(special) == 1):
        if len(special) == 1:
            return cpp_data.id_from_level(<int> nuc, <double> level, <std_string> spc)
        else:
            return cpp_data.id_from_level(<int> nuc, <double> level)
    else:
        return nuc


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
    nuc : int or None
        state_id of metastable state or None if the metastable state
        can not be found.
    """
    rtn = cpp_data.metastable_id(<int> nuc, <int> level)
    if rtn < 0:
        rtn = None
    return rtn


def decay_half_life(from_nuc, to_nuc):
    """
    Returns the half life from ENSDF decay dataset parent data

    Parameters
    ----------
    from_nuc : int
        parent nuclide in state_id form
    to_nuc : int
        child nuclide in state_id form

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
        parent nuclide in state_id form

    Returns
    -------
    half_lives : array of pairs
        An array of half lives and their errors for a given parent nuclide
    """
    half_lives = cpp_data.decay_half_lifes(<int> parent)
    return half_lives


def decay_data_children(parent):
    """
    Returns a list of child nuclides from the given parent based on decay
    datasets with the input parent listed.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    children : vector of ints
        An vector of child state_id's for a given parent nuclide
    """
    return cpp_data.decay_data_children(<int> parent)


def decay_branch_ratio(from_nuc, to_nuc):
    """
    Returns the branch ratio from ENSDF decay dataset data

    Parameters
    ----------
    from_nuc : int
        parent nuclide in state_id form
    to_nuc : int
        child nuclide in state_id form

    Returns
    -------
    ratio : double
        branching ratio
    """
    ratio, error = cpp_data.decay_branch_ratio(cpp_pair[int,int](from_nuc, to_nuc))
    return ratio, error

def decay_branch_ratio_byparent(parent):
    """
    Returns a list branching ratios from ENSDF decay dataset data for a given
    parent.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

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
        parent nuclide in state_id form
    to_nuc : int
        child nuclide in state_id form

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
        parent nuclide in state_id form

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
        parent nuclide in state_id form
    to_nuc : int
        child nuclide in state_id form

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
        parent nuclide in state_id form

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
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray energies and errors
    """
    return cpp_data.gamma_energy(<int> parent)


def gamma_energy_byen(en, enerror=None):
    """
    Returns a list of gamma ray energies from ENSDF decay dataset in a given
    energy range.

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
        An array of gamma ray energies and errors
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.gamma_energy(<double> en, <double> enerror)


def gamma_photon_intensity(parent):
    """
    Returns a list of gamma ray photon intensities from ENSDF decay dataset
    from a given parent

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray photon intensities and errors
    """
    return cpp_data.gamma_photon_intensity(<int> parent)


def gamma_photon_intensity_byen(en, enerror=None):
    """
    Returns a list of gamma ray photon intensities from ENSDF decay dataset
    from a given gamma ray energy

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
        An array of gamma ray photon intensities and errors
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.gamma_photon_intensity(<double> en,<double> enerror)

def gamma_conversion_intensity(parent):
    """
    Returns a list of gamma ray conversion intensities from ENSDF decay dataset
    from a given parent

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

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
        parent nuclide in state_id form

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
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of pairs
        An array of gamma ray level pairs in state_id form
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
        An array of gamma ray level pairs in state_id form
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.gamma_from_to(<double> en,<double> enerror)

def gamma_parent_child(en, enerror=None):
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
    ratios : array of pairs
        An array of gamma ray parents in state_id form
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.gamma_parent_child(<double> en, <double> enerror)


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
        An array of gamma ray parents in state_id form
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.gamma_parent(<double> en, <double> enerror)

def gamma_child_byen(en, enerror=None):
    """
    Returns a list of gamma ray children from ENSDF decay dataset
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
       An array of gamma ray children in state_id form
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.gamma_child(<double> en, <double> enerror)

def gamma_child_byparent(parent):
    """
    Returns a list of gamma ray children from ENSDF decay dataset
    based on gamma-ray energy.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form



    Returns
    -------
    ratios : array of ints
        An array of gamma ray children in state_id form
    """
    return cpp_data.gamma_child(<int> parent)

def gamma_xrays(parent):
    """
    Returns an array of arrays of xrays associated with the gamma
    rays from an input parent nuclide

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of arrays
        This returns an array of length 4 arrays containing pairs of energies
        and intensities of the following X-rays: Ka1, Ka2, Kb, L

    """
    return cpp_data.gamma_xrays(<int> parent)

def alpha_energy(parent):
    """
    Returns a list of alpha energies from ENSDF decay dataset from a given
    parent

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of alpha energies
    """
    return cpp_data.alpha_energy(<int> parent)

def alpha_intensity(parent):
    """
    Returns a list of alpha intensities from ENSDF decay dataset from a given
    parent

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of alpha intensities
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
        An array of alpha parents in state_id form
    """
    if enerror is None:
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
        An array of alpha children in state_id form
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.alpha_child(<double> en, <double> enerror)

def alpha_child_byparent(parent):
    """
    Returns a list of alpha children from ENSDF decay dataset
    based on alpha parent.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of ints
        An array of alpha children in state_id form
    """
    return cpp_data.alpha_child(<int> parent)

def beta_endpoint_energy(parent):
    """
    Returns a list of beta endpoint energies from ENSDF decay dataset
    based on parent nuclide.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of beta endpoint energies
    """
    return cpp_data.beta_endpoint_energy(<int> parent)

def beta_average_energy(parent):
    """
    Returns a list of beta average energies from ENSDF decay dataset
    based on parent nuclide.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of beta average energies
    """
    return cpp_data.beta_average_energy(<int> parent)

def beta_intensity(parent):
    """
    Returns a list of beta intensities from ENSDF decay dataset
    based on parent nuclide.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of beta intensities
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
    if enerror is None:
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
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.beta_child(<double> en, <double> enerror)

def beta_child_byparent(parent):
    """
    Returns a list of beta minus children from ENSDF decay dataset
    based on parent.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of ints
        An array of beta- children in state_id form
    """
    return cpp_data.beta_child(<int> parent)

def ecbp_endpoint_energy(parent):
    """
    Returns a list of beta plus endpoint energies from ENSDF decay dataset from
    a given parent.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of beta plus endpoint energies
    """
    return cpp_data.ecbp_endpoint_energy(<int> parent)

def ecbp_average_energy(parent):
    """
    Returns a list of beta plus average energies from ENSDF decay dataset from
    a given parent.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of beta plus average energies
    """
    return cpp_data.ecbp_average_energy(<int> parent)

def ec_intensity(parent):
    """
    Returns a list of electron capture intensisities from ENSDF decay dataset
    from a given parent.

    Parameters
    ----------
    parent : int
        parent nuclide in state_id form

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
        parent nuclide in state_id form

    Returns
    -------
    ratios : array of doubles
        An array of beta plus intensities
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
    if enerror is None:
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
        An array of beta plus/electron capture children in state_id form
    """
    if enerror is None:
        enerror = en * 0.01
    return cpp_data.ecbp_child(<double> en, <double> enerror)

def ecbp_child_byparent(parent):
    """
    Returns a list of beta plus children from ENSDF decay dataset
    based on parent.

    Parameters
    ----------
    parent : int
        parent nuclide in nuc_id form in state_id form

    Returns
    -------
    ratios : array of ints
        An array of beta+ children in state_id form
    """
    return cpp_data.ecbp_child(<int> parent)

def ecbp_xrays(parent):
    """
    Returns an array of arrays of xrays associated with the electron capture
    and beta plus decays from an input parent nuclide

    Parameters
    ----------
    parent : int
        parent nuclide

    Returns
    -------
    ratios : array of arrays
        This returns an array of length 4 arrays containing pairs of energies
        and intensities of the following X-rays: Ka1, Ka2, Kb, L

    """
    return cpp_data.ecbp_xrays(<int> parent)
