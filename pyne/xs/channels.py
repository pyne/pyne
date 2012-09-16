import collections

import numpy as np
# Hide warnings from numpy
np.seterr(divide='ignore')

import scipy.integrate
import tables as tb

import pyne
import pyne.data
import pyne.xs.models
from pyne import nucname
from pyne.xs import cache
from pyne.xs.models import group_collapse
from pyne.material import Material


def _prep_cache(xs_cache, E_g=None, phi_g=None):
    """Ensures that certain values are in the cache safely."""
    if E_g is not None:
        xs_cache['E_g'] = E_g

    if phi_g is not None:
        xs_cache['phi_g'] = phi_g


def _atom_weight_channel(chanfunc, nucspec, *args, **kwargs):
    """Convolves a channel for several nuclides based on atomic weights."""
    xs_cache = kwargs['xs_cache'] if 'xs_cache' in kwargs else cache.xs_cache
    # convert to atom weights
    if isinstance(nucspec, Material):
        aws = nucspec.to_atom_frac()
    elif isinstance(nucspec, collections.Mapping):
        aws = nucspec
    elif isinstance(nucspec, collections.Sequence):
        aws = dict(nucspec)

    # tally the channels as we go
    weight_total = 0.0
    chan = np.zeros(len(xs_cache['E_g']) - 1, float)
    for nuc, weight in aws.items():
        weight_total += weight
        nuc_chan = chanfunc(nuc, *args, **kwargs)
        chan += weight * nuc_chan

    # re-normalize
    if weight_total != 1.0:
        chan = chan / weight_total

    return chan


def sigma_f(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron fission cross-section for a nuclide for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.  If any of these are None-valued, 
    values from the cache are used.  The energy groups and fluxes are normally 
    ordered from highest-to-lowest energy.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        fission cross-section.
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    sigma_f_g : ndarray 
        An array of the collapsed fission cross-section.

    Notes
    -----
    This always pulls the fission cross-section out of nuc_data library.    

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collections.Iterable) and not isinstance(nuc, basestring):
        return _atom_weight_channel(sigma_f, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.zzaaam(nuc)
    key = (nuc, 'f', temp)
    return xs_cache[key]


def sigma_s_gh(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron scattering cross-section kernel for a nuclide for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that g, h index G, 
    n indexes N, and G < N.  g is for the incident energy and h is for the exiting energy.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        scattering kernel.
    T : float
        Tempurature of the target material [kelvin].
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    sig_s_gh : ndarray 
        An array of the scattering kernel.

    Notes
    -----
    This pulls the scattering length out of nuc_data library.

    Warnings
    --------
    This function is currently a stub until the proper way to compute the 
    scattering kernel is determined.  This function is safe to use but the 
    results are trivial.  This function simply returns an array with the 
    diagonal elements set to sigma_s as computed by pyne.xs.models.sigma_s().
    This conserves the calculation of sigma_s_g by summing sigma_s_gh over 
    the h-index.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collections.Iterable) and not isinstance(nuc, basestring):
        return _atom_weight_channel(sigma_s_gh, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.zzaaam(nuc)
    key = (nuc, 's_gh', temp)

    # Don't recalculate anything if you don't have to
    if key in xs_cache:
        return xs_cache[key]

    # Get some needed data
    E_g = xs_cache['E_g']
    G = len(E_g) - 1
    b = pyne.data.b(nuc)
    aw = pyne.data.atomic_mass(nuc)

    # OMG FIXME So hard!
    ## Initialize the scattering kernel
    #sig_s_gh = np.zeros((G, G), dtype=float)
    #
    ## Calculate all values of the kernel
    #for g, h in product(range(G), range(G)):
    #    # Numerator inetgration term 
    #    dnumer = lambda _E_prime, _E: sigma_s_E(_E, b, M_A, T) *  P(_E, _E_prime, M_A, T) * xs_cache['phi_g'][g]
    #
    #    # Integral
    #    nE = 26
    #    E_space = np.logspace(np.log10(xs_cache['E_g'][g]), np.log10(xs_cache['E_g'][g+1]), nE)
    #    E_prime_space = np.logspace(np.log10(xs_cache['E_g'][h]), np.log10(xs_cache['E_g'][h+1]), nE)
    #
    #    numer = msmintegrate.dbltrapz(dnumer, E_space, E_prime_space)
    #
    #    # Denominator term, analytically integrated
    #    denom = xs_cache['phi_g'][g] * (xs_cache['E_g'][g+1] - xs_cache['E_g'][g])
    #
    #    # Cross section value
    #    sig_s_gh[g, h] = numer / denom

    # Temporary stub
    E_g_centers = (E_g[1:] + E_g[:-1]) / 2.0
    sig_s = pyne.xs.models.sigma_s(E_g_centers, b, aw, temp)
    sig_s_gh = np.diag(sig_s)

    xs_cache[key] = sig_s_gh
    return sig_s_gh


def sigma_s(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron scattering cross-section for a nuclide. 

    .. math::
        \\sigma_{s, g} = \\sum_{h} \\sigma_{s, g\\to h} 

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        scattering cross section.
    T : float
        Tempurature of the target material [kelvin].
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    sig_s_g : ndarray 
        An array of the scattering cross section.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collections.Iterable) and not isinstance(nuc, basestring):
        return _atom_weight_channel(sigma_s, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.zzaaam(nuc)
    key_g  = (nuc, 's_g', temp)
    key_gh = (nuc, 's_gh', temp)

    # Don't recalculate anything if you don't have to
    if key_g in xs_cache:
        return xs_cache[key_g]

    # This calculation requires the scattering kernel
    if key_gh not in xs_cache:
        xs_cache[key_gh] = sigma_s_gh(nuc, temp, group_struct, phi_g, xs_cache)

    # Sum over all h
    sig_s_g = xs_cache[key_gh].sum(axis=1)

    # Put this value back into the cache, with the appropriate label
    xs_cache[key_g] = sig_s_g
    return sig_s_g


def sigma_a_reaction(nuc, rx, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron absorption reaction cross-section for a nuclide for a 
    new, lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        absorption reaction cross-section.
    rx : str
        Reaction key. ('gamma', 'alpha', 'p', etc.)
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    sigma_rx_g : ndarray 
        An array of the collapsed absorption reaction cross section.

    Notes
    -----
    This always pulls the absorption reaction cross-section out of the nuc_data.    

    See Also
    --------
    pyne.xs.cache.ABSORPTION_RX
    pyne.xs.cache.ABSORPTION_RX_MAP 
    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collections.Iterable) and not isinstance(nuc, basestring):
        return _atom_weight_channel(sigma_a_reaction, nuc, rx=rx, temp=temp, 
                                    xs_cache=xs_cache)
    nuc = nucname.zzaaam(nuc)
    key= (nuc, rx, temp)
    return xs_cache[key]


def metastable_ratio(nuc, rx, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the ratio between a reaction that leaves the nuclide in a 
    metastable state and the equivalent reaction that leaves the nuclide in 
    the ground state.  This allows the calculation of metastable cross sections 
    via sigma_ms = ratio * sigma_ground.  Note that g indexes G, n indexes N, 
    and G < N.

    Note: This always pulls the absorption reaction cross-sections out of the nuc_data library.    

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        metastable ratio.
    rx : str
        Reaction key. ('gamma', 'alpha', 'p', etc.)
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    ratio_rx_g : ndarray
        An array of the ratio of the metastable cross section for a reaction 
        to the ground state reaction.

    Notes
    -----
    This always pulls the absorption reaction cross section out of the nuc_data.

    See Also
    --------
    pyne.xs.cache.ABSORPTION_RX
    pyne.xs.cache.ABSORPTION_RX_MAP 
    """
    if isinstance(nuc, int) or isinstance(nuc, basestring):
        xs_cache = cache.xs_cache if xs_cache is None else xs_cache
        _prep_cache(xs_cache, group_struct, phi_g)
        nuc = nucname.zzaaam(nuc)
        key = (nuc, rx + '_x_ratio', temp)
        if key in xs_cache:
            return xs_cache[key]

    # Get the cross-sections
    sigma_rx = sigma_a_reaction(nuc, rx, temp, group_struct, phi_g, xs_cache)
    sigma_rx_x = sigma_a_reaction(nuc, rx + '_x', temp, group_struct, phi_g, xs_cache)

    # Get the ratio
    ratio_rx_g = sigma_rx_x / sigma_rx
    ratio_rx_g[ratio_rx_g < 0.0] = 0.0
    ratio_rx_g[ratio_rx_g == np.inf] = 0.0
    ratio_rx_g[np.isnan(ratio_rx_g)] = 0.0

    if isinstance(nuc, int):
        xs_cache[key] = ratio_rx_g
    return ratio_rx_g


def sigma_a(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron absorption cross section for a nuclide for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        absorption cross section.
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    sigma_a_g : ndarray 
        An array of the collapsed absorption cross section.

    Notes
    -----
    This always pulls the absorption cross section out of the nuc_data.    

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collections.Iterable) and not isinstance(nuc, basestring):
        return _atom_weight_channel(sigma_a, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.zzaaam(nuc)
    key = (nuc, 'a', temp)
    return xs_cache[key]


def chi(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None, eres=101):
    """Calculates the neutron fission energy spectrum for an isotope for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        neutron fission energy spectrum.
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.
    eres : int
        Number of energy-points to integrate over per group.

    Returns
    -------
    chi_g : ndarray 
        An array of the fission energy spectrum.

    See Also
    --------
    pyne.xs.models.chi : used under the covers by this function.
    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collections.Iterable) and not isinstance(nuc, basestring):
        return _atom_weight_channel(chi, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.zzaaam(nuc)
    key = (nuc, 'chi', temp)

    # Don't recalculate anything if you don't have to
    if key in xs_cache:
        return xs_cache[key]

    # Get the the set of nuclides we know we need chi for.  
    if 'fissionable_nucs' not in xs_cache:
        with tb.openFile(pyne.nuc_data, 'r') as f:
            if '/neutron/cinder_xs/fission' in f:
                fn = set(f.root.neutron.cinder_xs.fission.cols.nuc)
            else:
                fn = set()
        xs_cache['fissionable_nucs'] = fn
    fissionable_nucs = xs_cache['fissionable_nucs']
    if (nuc not in fissionable_nucs) and (86 <= nuc/10000):
        fissionable_nucs.add(nuc)

    # Perform the group collapse on a continuous chi
    E_g = xs_cache['E_g']
    G = len(E_g) - 1
    chi_g = np.zeros(G, dtype='f8')
    if (nuc in fissionable_nucs):
        for g in range(G):
            E_space = np.logspace(np.log10(E_g[g]), np.log10(E_g[g+1]), eres)
            dnumer = pyne.xs.models.chi(E_space)
            numer = scipy.integrate.trapz(dnumer, E_space)
            denom = (E_g[g+1] - E_g[g])
            chi_g[g] = (numer / denom)
        # renormalize chi
        chi_g = chi_g / chi_g.sum()
    # Put this value back into the cache, with the appropriate label
    xs_cache[key] = chi_g
    return chi_g



def sigma_t(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the total neutron cross section for a nuclide. 

    .. math::
        \\sigma_{t, g} = \\sigma_{a, g} + \\sigma_{s, g}

    Parameters
    ----------
    nuc : int, str, Material, or dict-like 
        A nuclide or nuclide-atom fraction mapping for which to calculate the 
        total cross section.
    T : float, optional
        Tempurature of the target material [kelvin].
    E_g : array-like of floats, optional
        New, lower fidelity energy group structure [MeV] that is of length G+1. 
    E_n : array-like of floats, optional
        Higher resolution energy group structure [MeV] that is of length N+1. 
    phi_n : array-like of floats, optional
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over.  
        Length N.  

    Returns
    -------
    sig_t_g : ndarray 
        An array of the total cross section.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collections.Iterable) and not isinstance(nuc, basestring):
        return _atom_weight_channel(sigma_t, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.zzaaam(nuc)
    key_a = (nuc, 'a', temp)
    key_s = (nuc, 's', temp)
    key_t = (nuc, 't', temp)

    # Don't recalculate anything if you don't have to
    if key_t in xs_cache:
        return xs_cache[key_t]

    # This calculation requires the abosorption cross-section
    if key_a not in xs_cache:
        xs_cache[key_a] = sigma_a(nuc, temp, group_struct, phi_g, xs_cache)
    if key_s not in xs_cache:
        xs_cache[key_s] = sigma_s(nuc, temp, group_struct, phi_g, xs_cache)
    sig_t_g = xs_cache[key_a] + xs_cache[key_s]
    xs_cache[key_t] = sig_t_g
    return sig_t_g


