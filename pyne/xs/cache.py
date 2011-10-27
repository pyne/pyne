"""This module provides a cross-section cache which automatically extracts cross-sections
from the nuclear database."""
from itertools import product
from collections import MutableMapping

import numpy as np
import tables as tb

from pyne import nucname
from pyne.pyne_config import pyne_conf


###############################################################################
### Set up a cross-section cache so the same data isn't loaded repetitively ###
###############################################################################

def is_g_indexed(key):
    if isinstance(key, basestring):
        is_g = '_g' in key
    else:
        is_g = '_g' in key[0]
    return is_g

class XSCache(MuatbleMapping):
    """A lightweight multigroup cross-section cache based off of python dictionaries.
    High resolution (*_n) data will be read from nuc_data.  Note, that this requires
    that nuc_data.h5 was built with CINDER data.
    """

    def __init__(self):
        self._cache = {}

        self._get_fns = {'E_n': get_E_n,
                         'sigma_f_n': get_sigma_f_n,
                         'sigma_a_n': get_sigma_a_n,
                         'sigma_rx_n': get_sigma_reaction_n,
                         'phi_g': lambda: get_phi_g(self['E_g'], self['E_n'], self['phi_n']),
                        }

    #
    # Mutable mapping pass-through interface
    #

    def __len__(self):
        return len(self._cache)

    def __iter__(self):
        return iter(self._cache)

    def __contains__(self, key):
        return (key in self._cache)

    def __delitem__(self, key):
        del self._cache[key]

    #
    # Explicit overrides
    #

    def __getitem__(self, key):
        """Key lookup by via custom loading from the nuc_data database file."""

        if key not in self._cache:
            if isinstance(key, basestring):
                self._cache[key] = self._get_fns[key]()
            else:
                self._cache[key] = self._get_fns[key[0]](*key[1:])

        # Return the value requested
        return self._cache[key]


    def __setitem__(self, key, value):
        """Key setting via custom cache functionality."""
        # Set the E_g
        if (key == 'E_g'):
            value = np.array(value, dtype=float)

            # If the E_gs are the same, don't set anything
            if ('E_g' in self):
                if (len(value) == len(self['E_g'])) and (value == self['E_g']).all():
                    return 

            # Otherwise, preload some stuff.
            self._cache['partial_energy_matrix'] = partial_energy_matrix(value, self['E_n'])

            # And remove any previous paramters dependent on E_g
            dirty_keys = [k for k in self._cache if is_g_indexed(k)]
            for dk in dirty_keys:
                del self._cache[dk]


        # Set the E_n
        if (key == 'E_n'):
            value = np.array(value, dtype=float)
            
            # If the E_gs are the same, don't set anything
            if ('E_n' in self):
                if (len(value) == len(self['E_n'])) and (value == self['E_n']).all():
                    return 

            # Otherwise, preload some stuff.
            if 'E_g' in self:
                self._cache['partial_energy_matrix'] = partial_energy_matrix(self['E_g'], value)


        # Set the high resolution flux, phi_n
        if key == 'phi_n':
            value = np.array(value, dtype=float)

            # If the flux is same, don't set or remove anything
            if ('phi_n' in self):
                if (len(value) == len(self['phi_n'])) and (value == self['phi_n']).all():
                    return 

            # And remove any previous paramters dependent on phi_n
            dirty_keys = [k for k in self._cache if (is_g_indexed(k) and k != 'E_g')]
            for dk in dirty_keys:
                del self._cache[dk]

        # Set the value normally
        self._cache[key] = value


#################################
### Get functions for XSCache ### 
#################################

def get_E_n():
    """Loads and returns the cinder energy bounds array, E_n, from the nuc_data library."""
    nuc_data = pyne_conf.NUC_DATA_PATH
    with tb.openFile(nuc_data, 'r') as f:
        E_n = np.array(f.root.neutron.cinder_xs.E_g)
    return E_n


def get_sigma_f_n(nuc):
    """Grabs a nuclide's fission cross-section from the nuc_data library.

    Parameters
    ----------
    nuc : int
        nuclide in zzaaam form.

    Returns
    -------
    sigma_f_n : numpy array 
        This nuclide's fission cross-section pulled from the nuc_data library file.  If not 
        present in the library, a zero-array is returned.
    """
    nuc_data = pyne_conf.NUC_DATA_PATH
    with tb.openFile(nuc_data, 'r') as f:
        N = f.root.neutron.cinder_xs.fission.coldescrs['xs'].shape[0]
        cond = 'nuc_zz == {0}'.format(nuc)
        rows = [np.array(row['xs']) for row in f.root.neutron.cinder_xs.fission.where(cond)]

    if len(rows) == 0:
        # Not fissionable, return zero-array
        sigma_f_n = np.zeros(N, dtype=float)
    elif len(rows) == 1:
        # Return cross-section from file
        sigma_f_n = rows[0]
    else:
        nuc_name = nucname.name(nuc)
        err_str = "The database contains multiple entries for the fission cross-section for {0}!".format(nuc_name)
        raise ValueError(err_str)

    return sigma_f_n


def get_sigma_a_n(nuc):
    """Grabs a nuclide's absorption cross-section from the nuc_data library.

    Parameters
    ----------
    nuc : int
        nuclide in zzaaam form.

    Returns
    -------
    sigma_a_n : numpy array
        This nuclide's absorption cross-section pulled from the database library file.  If not 
        present in the library, a zero-array is returned.
    """
    nuc_data = pyne_conf.NUC_DATA_PATH
    with tb.openFile(nuc_data, 'r') as f:
        N = f.root.neutron.cinder_xs.absorption.coldescrs['xs'].shape[0]
        cond = "(from_nuc_zz == {0}) & (reaction_type != 'c')".format(nuc)
        rows = [np.array(row['xs']) for row in f.root.neutron.cinder_xs.absorption.where(cond)]

    if len(rows) == 0:
        # No absportion, return zero-array
        sigma_a_n = np.zeros(N, dtype=float)
    else:
        rows = np.array(rows)
        sigma_a_n = rows.sum(axis=0)

    # Add in the fission cross-section
    sigma_f_n = get_sigma_f_n(nuc)
    sigma_a_n += sigma_f_n

    return sigma_a_n



ABSORPTION_RX = set(["", 'np *', 'a  *', 'h  *', '2p *', '3n *', 'd  *', 'np/d', 
                     'na', '*', 'nd', 'g  *', '3n', 'np', 'nt', 't', 'nt *', 
                     'x  *', '4n *', 'na *', 'nd *', 't  *', 'a', 'c', '2p', 'd', 
                     'g', 'h', 'n', '4n', 'p', 'n  *', '2a', '2n *', 'x', '2n', 
                     'nh *', 'p  *'])

ABSORPTION_RX_MAP = {
    'neutron': 'n',
    'gamma': 'g', 
    'alpha': 'a',
    'proton': 'p',
    'trit': 't',
    'triton': 't',
    'deut': 'd',
    'deuteron': 'd',
    }

def get_sigma_reaction_n(nuc, rx):
    """Grabs a nuclide's absorption reaction cross-section from the nuc_data library.

    Parameters
    ----------
    nuc : int
        nuclide in zzaaam form.
    rx : str 
        Reaction key: 'gamma', 'alpha', 'p', etc.

    Returns
    -------
    sigma_rx_n : numpy array
        This nuclide's absorption reaction cross-section pulled from the nuc_data library 
        file.  If not present in the library, a zero-array is returned.
    """
    # munge the reaction rate
    if rx not in ABSORPTION_RX:
        while any([(key in rx) for key in ABSORPTION_RX_MAP]):
            for key, value in ABSORPTION_RX_MAP.items():
                rx = rx.replace(key, value)

    if '_x' in rx:
        if len(rx) == 3:
            rx = rx.replace('_x', '  *')
        else:
            rx = rx.replace('_x', ' *')

    if rx not in ABSORPTION_RX:
        msg = "the reaction '{rx}' is not valid.".format(rx=rx)
        raise IndexError(msg)

    # Read in the data
    nuc_data = pyne_conf.NUC_DATA_PATH
    with tb.openFile(nuc_data, 'r') as f:
        N = f.root.neutron.cinder_xs.absorption.coldescrs['xs'].shape[0]
        cond = "(from_nuc_zz == {0}) & (reaction_type == '{1}')".format(nuc, rx)
        rows = [np.array(row['xs']) for row in f.root.neutron.cinder_xs.absorption.where(cond)]

    if len(rows) == 0:
        # No absportion, return zero-array
        sigma_rx_n = np.zeros(N, dtype=float)
    else:
        rows = np.array(rows)
        sigma_rx_n = rows.sum(axis=0)

    return sigma_rx_n



##############################
### Partial group collapse ###
##############################

def partial_energy_matrix(E_g, E_n):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.

    Args:
        * E_n (sequence of floats): higher resolution energy group structure [MeV]
          that is of length N+1. Ordered from lowest-to-highest energy.
        * E_g (sequence of floats): lower resolution energy group structure [MeV]
          that is of length G+1. Ordered from lowest-to-highest energy.

    Returns:
        * pem (array of fractions): This is a GxN sized matrix that when dotted with a 
          high-resolution flux produces a low-resolution flux.
    """
    # Some convienence paramters
    G = len(E_g) - 1
    N = len(E_n) - 1

    index_E_n = np.arange(N+1)

    # Some assertions to ensure that everything is well formed
    assert E_n[0] <= E_g[0]
    assert E_g[-1] <= E_n[-1]

    # Get the interior points for each gth group in n-space
    inner_mask = np.array([(E_g[g] <= E_n) & (E_n <= E_g[g+1]) for g in range(G)])

    # Get the upper and lower nth index for every gth group
    lower_index = np.array([index_E_n[inner_mask[g]][0] for g in range(G)])
    upper_index = np.array([index_E_n[inner_mask[g]][-1] for g in range(G)])

    # Convert the mask to initialize the partial enery matrix
    # Hack off the last index of the mask to make the right size
    pem = np.array(inner_mask[:, :-1], dtype=float)

    # Check for partial contibutions at the edges
    for g in range(G):
        # Lower bound
        lig = lower_index[g]
        if lig != 0:
            pem[g][lig-1] = (E_n[lig] - E_g[g]) / (E_n[lig] - E_n[lig-1])

        # Upper bound
        uig = upper_index[g]
        if uig < N:
            pem[g][uig] = (E_g[g+1] - E_n[uig]) / (E_n[uig+1] - E_n[uig])

    return pem


def get_phi_g(E_g, E_n, phi_n):
    """Calculates the lower resolution flux, phi_g, from the lower resolution group stucture E_g, 
    the higher resolution groups E_n, and the higher resolution flux phi_n.

    Parameters
    ----------
    E_g : sequence of floats 
        Lower resolution energy group structure [MeV] that is of length G+1. 
        Ordered from lowest-to-highest energy.
    E_n : sequence of floats 
        Higher resolution energy group structure [MeV] that is of length N+1. 
        Ordered from lowest-to-highest energy.
    phi_n : sequence of floats
        The high-fidelity flux [n/cm^2/s] to collapse the fission cross-section over (length N).  
        Ordered from lowest-to-highest energy.

    Returns
    -------
    phi_g : numpy array of floats 
        The flux collapsed to G energy groups.
    """
    pem = partial_energy_matrix(E_g, E_n)
    phi_g = np.dot(pem, phi_n)
    return phi_g


xs_cache = XSCache()
