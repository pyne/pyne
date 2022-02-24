"""This module provides an easy interface to very quickly grab multigroup cross 
sections from the cross section cache and collapse them to the appropriate group 
structure.  Additionally, it provides interfaces for some higher level functionality, 
such as computing cross sections for materials, fission energy spectra, metastable
ratios, etc.
"""
from __future__ import division
import sys
import collections

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections
from pyne.utils import QA_warn

import numpy as np
import scipy.integrate
import tables as tb

from .. import nuc_data
from .. import data
from .. import nucname
from .. import rxname
from ..material import Material
from . import models
from . import cache
from .models import group_collapse

QA_warn(__name__)

if sys.version_info[0] > 2:
    basestring = str

np.seterr(all="ignore")


def _prep_cache(xs_cache, E_g=None, phi_g=None):
    """Ensures that certain values are in the cache safely."""
    if E_g is not None:
        xs_cache["E_g"] = E_g

    if phi_g is not None:
        xs_cache["phi_g"] = phi_g


def _atom_mass_channel(chanfunc, nucspec, *args, **kwargs):
    """Convolves a channel for several nuclides based on atomic mass."""
    xs_cache = kwargs["xs_cache"] if "xs_cache" in kwargs else cache.xs_cache
    # convert to atomic mass
    if isinstance(nucspec, Material):
        aws = nucspec.to_atom_frac()
    elif isinstance(nucspec, collectionsAbc.Mapping):
        aws = nucspec
    elif isinstance(nucspec, collectionsAbc.Sequence):
        aws = dict(nucspec)

    # tally the channels as we go
    mass_total = 0.0
    chan = np.zeros(len(xs_cache["E_g"]) - 1, float)
    for nuc, mass in aws.items():
        mass_total += mass
        nuc_chan = chanfunc(nuc, *args, **kwargs)
        chan += mass * nuc_chan

    # re-normalize
    if mass_total != 1.0:
        chan = chan / mass_total

    return chan


def sigma_f(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron fission cross section for a nuclide.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like
        A nuclide or nuclide-atom fraction mapping.
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.

    Returns
    -------
    sigma_f_g : ndarray
        The fission cross-section, length G.

    Notes
    -----
    This always pulls the fission cross section out of the cache.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collectionsAbc.Iterable) and not isinstance(nuc, basestring):
        return _atom_mass_channel(sigma_f, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.id(nuc)
    key = (nuc, rxname.id("fission"), temp)
    return xs_cache[key]


def sigma_s_gh(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron scattering kernel for a nuclide.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like
        A nuclide or nuclide-atom fraction mapping.
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.

    Returns
    -------
    sig_s_gh : ndarray
        The scattering kernel, shape GxG.

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
    if isinstance(nuc, collectionsAbc.Iterable) and not isinstance(nuc, basestring):
        return _atom_mass_channel(sigma_s_gh, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.id(nuc)
    key = (nuc, "s_gh", temp)

    # Don't recalculate anything if you don't have to
    if key in xs_cache:
        return xs_cache[key]

    # Get some needed data
    E_g = xs_cache["E_g"]
    G = len(E_g) - 1
    b = data.b(nuc)
    aw = data.atomic_mass(nuc)

    # OMG FIXME So hard!
    ## Initialize the scattering kernel
    # sig_s_gh = np.zeros((G, G), dtype=float)
    #
    ## Calculate all values of the kernel
    # for g, h in product(range(G), range(G)):
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
    sig_s = models.sigma_s(E_g_centers, b, aw, temp)
    sig_s_gh = np.diag(sig_s)

    xs_cache[key] = sig_s_gh
    return sig_s_gh


def sigma_s(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron scattering cross section for a nuclide.

    .. math::
        \\sigma_{s, g} = \\sum_{h} \\sigma_{s, g\\to h}

    Parameters
    ----------
    nuc : int, str, Material, or dict-like
        A nuclide or nuclide-atom fraction mapping.
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.

    Returns
    -------
    sig_s_g : ndarray
        The scattering cross section, length G.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collectionsAbc.Iterable) and not isinstance(nuc, basestring):
        return _atom_mass_channel(sigma_s, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.id(nuc)
    key_g = (nuc, "s_g", temp)
    key_gh = (nuc, "s_gh", temp)

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
    """Calculates the neutron absorption reaction cross section for a nuclide.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like
        A nuclide or nuclide-atom fraction mapping.
    rx : str
        Reaction key. ('gamma', 'alpha', 'p', etc.)
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.

    Returns
    -------
    sigma_rx_g : ndarray
        The collapsed absorption reaction cross section, length G.

    Notes
    -----
    This always pulls the absorption reaction cross-section out of the nuc_data.

    See Also
    --------
    pyne.xs.data_source.RX_TYPES
    pyne.xs.data_source.RX_TYPES_MAP

    """
    rx = rxname.id(rx)
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collectionsAbc.Iterable) and not isinstance(nuc, basestring):
        return _atom_mass_channel(
            sigma_a_reaction, nuc, rx=rx, temp=temp, xs_cache=xs_cache
        )
    nuc = nucname.id(nuc)
    key = (nuc, rx, temp)
    return xs_cache[key]


def metastable_ratio(nuc, rx, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the ratio between a reaction that leaves the nuclide in a
    metastable state and the equivalent reaction that leaves the nuclide in
    the ground state.  This allows the calculation of metastable cross sections
    via sigma_ms = ratio * sigma_ground.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like
        A nuclide or nuclide-atom fraction mapping.
    rx : str
        Reaction key. ('gamma', 'alpha', 'p', etc.)
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.

    Returns
    -------
    ratio_rx_g : ndarray
        An array of the ratio of the metastable cross section for a reaction
        to the ground state reaction.

    Notes
    -----
    This always pulls the absorption reaction cross section out of the cache.

    See Also
    --------
    pyne.xs.data_source.RX_TYPES
    pyne.xs.data_source.RX_TYPES_MAP

    """
    if isinstance(nuc, int) or isinstance(nuc, basestring):
        xs_cache = cache.xs_cache if xs_cache is None else xs_cache
        _prep_cache(xs_cache, group_struct, phi_g)
        nuc = nucname.id(nuc)
        key = (nuc, rx + "_x_ratio", temp)
        if key in xs_cache:
            return xs_cache[key]

    # Get the cross-sections
    sigma_rx = sigma_a_reaction(nuc, rx, temp, group_struct, phi_g, xs_cache)
    sigma_rx_x = sigma_a_reaction(nuc, rx + "_1", temp, group_struct, phi_g, xs_cache)

    # Get the ratio
    ratio_rx_g = sigma_rx_x / sigma_rx
    ratio_rx_g[ratio_rx_g < 0.0] = 0.0
    ratio_rx_g[ratio_rx_g == np.inf] = 0.0
    ratio_rx_g[np.isnan(ratio_rx_g)] = 0.0

    if isinstance(nuc, int):
        xs_cache[key] = ratio_rx_g
    return ratio_rx_g


def sigma_a(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None):
    """Calculates the neutron absorption cross section for a nuclide.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like
        A nuclide or nuclide-atom fraction mapping.
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.

    Returns
    -------
    sigma_a_g : ndarray
        The collapsed absorption cross section, length G.

    Notes
    -----
    This always pulls the absorption cross section out of the cache.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collectionsAbc.Iterable) and not isinstance(nuc, basestring):
        return _atom_mass_channel(sigma_a, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.id(nuc)
    key = (nuc, rxname.id("absorption"), temp)
    return xs_cache[key]


def chi(nuc, temp=300.0, group_struct=None, phi_g=None, xs_cache=None, eres=101):
    """Calculates the neutron fission energy spectrum for a nuclide.

    Parameters
    ----------
    nuc : int, str, Material, or dict-like
        A nuclide or nuclide-atom fraction mapping.
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.
    eres : int, optional
        Number of energy-points to integrate over per group.

    Returns
    -------
    chi_g : ndarray
        The fission energy spectrum, length G.

    See Also
    --------
    pyne.xs.models.chi : used under the covers by this function.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collectionsAbc.Iterable) and not isinstance(nuc, basestring):
        return _atom_mass_channel(chi, nuc, temp=temp, xs_cache=xs_cache, eres=eres)
    nuc = nucname.id(nuc)
    key = (nuc, "chi", temp)

    # Don't recalculate anything if you don't have to
    if key in xs_cache:
        return xs_cache[key]

    # Get the the set of nuclides we know we need chi for.
    if "fissionable_nucs" not in xs_cache:
        with tb.open_file(nuc_data, "r") as f:
            if "/neutron/cinder_xs/fission" in f:
                fn = set(f.root.neutron.cinder_xs.fission.cols.nuc)
            else:
                fn = set()
        xs_cache["fissionable_nucs"] = fn
    fissionable_nucs = xs_cache["fissionable_nucs"]
    if (nuc not in fissionable_nucs) and (86 <= nucname.znum(nuc)):
        fissionable_nucs.add(nuc)

    # Perform the group collapse on a continuous chi
    E_g = xs_cache["E_g"]
    G = len(E_g) - 1
    chi_g = np.zeros(G, dtype="f8")
    if nuc in fissionable_nucs:
        for g in range(G):
            E_space = np.logspace(np.log10(E_g[g]), np.log10(E_g[g + 1]), eres)
            dnumer = models.chi(E_space)
            numer = scipy.integrate.trapz(dnumer, E_space)
            denom = E_g[g + 1] - E_g[g]
            chi_g[g] = numer / denom
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
        A nuclide or nuclide-atom fraction mapping.
    temp : float, optional
        Temperature [K] of material, defaults to 300.0.
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1,
        defaults to xs_cache['E_g'].
    phi_g : array-like of floats, optional
        Group fluxes [n/cm^2/s] matching group_struct, length G, defaults to
        xs_cache['phi_g'].
    xs_cache : XSCache, optional
        Cross section cache to use, defaults to pyne.xs.cache.xs_cache.

    Returns
    -------
    sig_t_g : ndarray
        The total cross section, length G.

    """
    xs_cache = cache.xs_cache if xs_cache is None else xs_cache
    _prep_cache(xs_cache, group_struct, phi_g)
    if isinstance(nuc, collectionsAbc.Iterable) and not isinstance(nuc, basestring):
        return _atom_mass_channel(sigma_t, nuc, temp=temp, xs_cache=xs_cache)
    nuc = nucname.id(nuc)
    key_a = (nuc, rxname.id("absorption"), temp)
    key_s = (nuc, rxname.id("scattering"), temp)
    key_t = (nuc, rxname.id("total"), temp)

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
