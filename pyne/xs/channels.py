import numpy as np
import scipy.integrate
import tables as tb

import pyne.data
import pyne.xs.models
from pyne import nucname
from pyne.xs.cache import xs_cache
from pyne.xs.models import group_collapse
from pyne.pyne_config import pyne_conf

# Hide warnings from numpy
np.seterr(divide='ignore')


def _prep_cache(E_g=None, E_n=None, phi_n=None):
    """Ensures that certain values are in the cache safely."""
    if E_n is not None:
        xs_cache['E_n'] = E_n

    # needs to follow E_n to calc collapse matrix properly
    if E_g is not None:
        xs_cache['E_g'] = E_g

    if phi_n is not None:
        xs_cache['phi_n'] = phi_n


def sigma_f(nuc, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron fission cross-section for a nuclide for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.  If any of these are None-valued, 
    values from the cache are used.  The energy groups and fluxes are normally 
    ordered from highest-to-lowest energy.

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the fission cross-section.
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
    _prep_cache(E_g, E_n, phi_n)

    # Get the fission XS
    nuc_zz = nucname.zzaaam(nuc)
    sigma_f_n_nuc_zz = ('sigma_f_n', nuc_zz)
    sigma_f_g_nuc_zz = ('sigma_f_g', nuc_zz)

    # Don't recalculate anything if you don't have to
    if sigma_f_g_nuc_zz in xs_cache:
        return xs_cache[sigma_f_g_nuc_zz]
    else:
        sigma_f_n = xs_cache[sigma_f_n_nuc_zz]

    # Perform the group collapse, knowing that the right data is in the cache
    sigma_f_g = group_collapse(sigma_f_n, xs_cache['phi_n'], phi_g=xs_cache['phi_g'], 
                               partial_energies=xs_cache['partial_energy_matrix'])

    # Put this value back into the cache, with the appropriate label
    xs_cache[sigma_f_g_nuc_zz] = sigma_f_g

    return sigma_f_g


def sigma_s_gh(nuc, T, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron scattering cross-section kernel for a nuclide for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that g, h index G, 
    n indexes N, and G < N.  g is for the incident energy and h is for the exiting energy.

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the scattering kernel.
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
    _prep_cache(E_g, E_n, phi_n)

    nuc_zz = nucname.zzaaam(nuc)
    key = ('sigma_s_gh', nuc_zz, T)

    # Don't recalculate anything if you don't have to
    if key in xs_cache:
        return xs_cache[key]

    # Get some needed data
    G = len(xs_cache['E_g']) - 1
    b = pyne.data.b(nuc_zz)
    aw = pyne.data.nuc_weight(nuc_zz)

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
    E_g = xs_cache['E_g']
    E_g_centers = (E_g[1:] + E_g[:-1]) / 2.0
    sig_s = pyne.xs.models.sigma_s(E_g_centers, b, aw, T)
    sig_s_gh = np.diag(sig_s)

    # Put this value back into the cache, with the appropriate label
    xs_cache[key] = sig_s_gh

    return sig_s_gh


def sigma_s(nuc, T, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron scattering cross-section for a nuclide. 

    .. math::
        \\sigma_{s, g} = \\sum_{h} \\sigma_{s, g\\to h} 

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the scattering cross section.
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
    _prep_cache(E_g, E_n, phi_n)

    nuc_zz = nucname.zzaaam(nuc)
    key_g = ('sigma_s_g', nuc_zz, T)
    key_gh = ('sigma_s_gh', nuc_zz, T)

    # Don't recalculate anything if you don't have to
    if key_g in xs_cache:
        return xs_cache[key_g]

    # This calculation requires the scattering kernel
    if key_gh not in xs_cache:
        xs_cache[key_gh] = sigma_s_gh(nuc, T, E_g, E_n, phi_n)

    # Sum over all h
    sig_s_g = xs_cache[key_gh].sum(axis=1)

    # Put this value back into the cache, with the appropriate label
    xs_cache[key_g] = sig_s_g

    return sig_s_g


def sigma_a_reaction(nuc, rx, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron absorption reaction cross-section for a nuclide for a 
    new, lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the absorption reaction cross-section.
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
    _prep_cache(E_g, E_n, phi_n)

    # Get the absorption XS
    nuc_zz = nucname.zzaaam(nuc)
    key_n = ('sigma_rx_n', nuc_zz, rx)
    key_g = ('sigma_rx_g', nuc_zz, rx)

    # Don't recalculate anything if you don't have to
    if key_g in xs_cache:
        return xs_cache[key_g]
    else:
        sigma_rx_n = xs_cache[key_n]

    # Perform the group collapse, knowing that the right data is in the cache
    sigma_rx_g = group_collapse(sigma_rx_n, xs_cache['phi_n'], phi_g=xs_cache['phi_g'], 
                                partial_energies=xs_cache['partial_energy_matrix'])

    # Put this value back into the cache, with the appropriate label
    xs_cache[key_g] = sigma_rx_g

    return sigma_rx_g


def metastable_ratio(nuc, rx, E_g=None, E_n=None, phi_n=None):
    """Calculates the ratio between a reaction that leaves the nuclide in a 
    metastable state and the equivalent reaction that leaves the nuclide in 
    the ground state.  This allows the calculation of metastable cross sections 
    via sigma_ms = ratio * sigma_ground.  Note that g indexes G, n indexes N, 
    and G < N.

    Note: This always pulls the absorption reaction cross-sections out of the nuc_data library.    

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the metastable.
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
    # Get the cross-sections
    sigma_rx = sigma_a_reaction(nuc, rx, E_g, E_n, phi_n)
    sigma_rx_x = sigma_a_reaction(nuc, rx + '_x', E_g, E_n, phi_n)

    # Get the ratio
    ratio_rx_g = sigma_rx_x / sigma_rx
    ratio_rx_g[ratio_rx_g < 0.0] = 0.0
    ratio_rx_g[ratio_rx_g == np.inf] = 0.0
    ratio_rx_g[np.isnan(ratio_rx_g)] = 0.0

    return ratio_rx_g


def sigma_a(nuc, E_g=None, E_n=None, phi_n=None):
    """Calculates the neutron absorption cross section for a nuclide for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the absorption cross section.
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
    _prep_cache(E_g, E_n, phi_n)

    # Get the absorption XS
    nuc_zz = nucname.zzaaam(nuc)
    key_n = ('sigma_a_n', nuc_zz)
    key_g = ('sigma_a_g', nuc_zz)

    # Don't recalculate anything if you don't have to
    if key_g in xs_cache:
        return xs_cache[key_g]
    else:
        sigma_a_n = xs_cache[key_n]

    # Perform the group collapse, knowing that the right data is in the cache
    sigma_a_g = group_collapse(sigma_a_n, xs_cache['phi_n'], phi_g=xs_cache['phi_g'], 
                               partial_energies=xs_cache['partial_energy_matrix'])

    # Put this value back into the cache, with the appropriate label
    xs_cache[key_g] = sigma_a_g

    return sigma_a_g


def chi(nuc, E_g=None, E_n=None, phi_n=None, eres=101):
    """Calculates the neutron fission energy spectrum for an isotope for a new, 
    lower resolution group structure using a higher fidelity flux.  Note that 
    g indexes G, n indexes N, and G < N.

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the neutron fission energy spectrum.
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
    _prep_cache(E_g, E_n, phi_n)

    # Get the fission XS
    nuc_zz = nucname.zzaaam(nuc)
    key = ('chi_g', nuc_zz)

    # Don't recalculate anything if you don't have to
    if key in xs_cache:
        return xs_cache[key]

    # Get the the set of nuclides we know we need chi for.  
    if 'fissionable_nucs' not in xs_cache:
        with tb.openFile(pyne_conf.NUC_DATA_PATH, 'r') as f:
            fn = set(f.root.neutron.cinder_xs.fission.cols.nuc_zz)
        xs_cache['fissionable_nucs'] = fn
    fissionable_nucs = xs_cache['fissionable_nucs']

    if (nuc_zz not in fissionable_nucs) and (86 <= nuc_zz/10000):
        fissionable_nucs.add(nuc_zz)

    # Perform the group collapse on a continuous chi
    G = len(xs_cache['E_g']) - 1
    chi_g = np.zeros(G, dtype=float)

    if (nuc_zz in fissionable_nucs):
        for g in range(G):
            E_space = np.logspace(np.log10(xs_cache['E_g'][g]), 
                                  np.log10(xs_cache['E_g'][g+1]), eres)
            dnumer = pyne.xs.models.chi(E_space)

            numer = scipy.integrate.trapz(dnumer, E_space)
            denom = (xs_cache['E_g'][g+1] - xs_cache['E_g'][g])

            chi_g[g] = (numer / denom)

        # renormalize chi
        chi_g = chi_g / chi_g.sum()

    # Put this value back into the cache, with the appropriate label
    xs_cache[key] = chi_g
    return chi_g



def sigma_t(nuc, T=300.0, E_g=None, E_n=None, phi_n=None):
    """Calculates the total neutron cross section for a nuclide. 

    .. math::
        \\sigma_{t, g} = \\sigma_{a, g} + \\sigma_{s, g}

    Parameters
    ----------
    nuc 
        A nuclide name for which to calculate the total cross section.
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
    _prep_cache(E_g, E_n, phi_n)

    # Get the total XS
    nuc_zz = nucname.zzaaam(nuc)
    key_a = ('sigma_a_g', nuc_zz)
    key_s = ('sigma_t_g', nuc_zz, T)
    key_t = ('sigma_t_g', nuc_zz, T)

    # Don't recalculate anything if you don't have to
    if key_t in xs_cache:
        return xs_cache[key_t]

    # This calculation requires the abosorption cross-section
    if key_a not in xs_cache:
        xs_cache[key_a] = sigma_a(nuc, E_g, E_n, phi_n)

    # This calculation requires the scattering cross-section
    if key_s not in xs_cache:
        xs_cache[key_s] = sigma_s(nuc, T, E_g, E_n, phi_n)

    # Sum over all h indeces
    sig_t_g = xs_cache[key_a] + xs_cache[key_s]

    # Put this value back into the cache, with the appropriate label
    xs_cache[key_t] = sig_t_g
    return sig_t_g


