"""This module provides a cross-section cache which automatically extracts cross-sections
from the nuclear database."""
from itertools import product
from collections import MutableMapping

import numpy as np
import tables as tb

from pyne import nucname
from pyne.pyne_config import pyne_conf
from pyne.xs.models import partial_energy_matrix, phi_g


###############################################################################
### Set up a cross-section cache so the same data isn't loaded repetitively ###
###############################################################################

def is_g_indexed(key):
    if isinstance(key, basestring):
        is_g = '_g' in key
    else:
        is_g = '_g' in key[0]
    return is_g

class XSCache(MutableMapping):
    """A lightweight multigroup cross-section cache based off of python dictionaries.
    High resolution (``*_n``) data will be read from nuc_data.  Note, that this requires
    that nuc_data.h5 was built with CINDER data.
    """

    def __init__(self):
        self._cache = {}

        self._get_fns = {'E_n': get_E_n,
                         'sigma_f_n': get_sigma_f_n,
                         'sigma_a_n': get_sigma_a_n,
                         'sigma_rx_n': get_sigma_a_reaction_n,
                         'phi_g': lambda: phi_g(self['E_g'], self['E_n'], self['phi_n']),
                         'eaf_xs': get_eaf_xs,
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

        if (key not in self._cache) and (key in self._get_fns or key[0] in self._get_fns):
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

def get_sigma_a_reaction_n(nuc, rx):
    """Grabs a nuclide's absorption reaction cross-section from the nuc_data library.

    Parameters
    ----------
    nuc : int
        Nuclide in zzaaam form.
    rx : str 
        Reaction key: 'gamma', 'alpha', 'p', etc.

    Returns
    -------
    sigma_rx_n : numpy array
        This nuclide's absorption reaction cross-section pulled from the nuc_data library 
        file.  If not present in the library, a zero-array is returned.

    See Also
    --------
    pyne.xs.cache.ABSORPTION_RX
    pyne.xs.cache.ABSORPTION_RX_MAP
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


#EAF_RX = set(["", 'np *', 'a  *', 'h  *', '2p *', '3n *', 'd  *', 'np/d', 
#                     'na', '*', 'nd', 'g  *', '3n', 'np', 'nt', 't', 'nt *', 
#                     'x  *', '4n *', 'na *', 'nd *', 't  *', 'a', 'c', '2p', 'd', 
#                     'g', 'h', 'n', '4n', 'p', 'n  *', '2a', '2n *', 'x', '2n', 
#                     'nh *', 'p  *',
EAF_RX = set(['1020', '1021', '1022', '1030', '1031', '1032', '1040', 
                     '1041', '1042', '1050', '1051', '1052', '1060', '1061', 
                     '1062', '1070', '1071', '1072', '1080', '1110', '1111', 
                     '1112', '160', '161', '162', '170', '171', '172', '180',
                     '220', '221', '222', '240', '280', '281', '282', '290',
                     '320', '321', '322', '330', '331', '332', '340', '341',
                     '342', '370', '371', '40', '41', '42', '420'])

EAF_RX_MAP = {
    'neutron': 'n',
    'gamma': 'g', 
    'alpha': 'a',
    'proton': 'p',
    'trit': 't',
    'triton': 't',
    'deut': 'd',
    'deuteron': 'd',
    }


def get_eaf_xs(nuc, rx):
    """Grabs group-wise cross section values fro a given nuclide and reaction from the EAF format activation library.

    Parameters
    ----------
    nuc : int
        Nuclide in zzaaam form.
    rx : str 
        Reaction MT # in nnnm form.
      OR: (eventually
        Reaction key: 'gamma', 'alpha', 'p', etc.

    Returns
    -------
    sigma_n : numpy array
        The cross-section pulled from the nuc_data library file.
        If not present in the library, a zero-array is returned.

    See Also
    --------
    pyne.xs.cache.EAF_RX
    pyne.xs.cache.EAF_RX_MAP
    """
    # munge the reaction rate
    #TODO:
    #if rx not in EAF_RX:
    #    while any([(key in rx) for key in EAF_RX_MAP]):
    #        for key, value in EAF_RX_MAP.items():
    #            rx = rx.replace(key, value)

    #if '_x' in rx:
    #    if len(rx) == 3:
    #        rx = rx.replace('_x', '  *')
    #    else:
    #        rx = rx.replace('_x', ' *')

    if rx not in EAF_RX:
        msg = "the reaction '{rx}' is not valid.".format(rx=rx)
        raise IndexError(msg)

    # Read in the data
    nuc_data = pyne_conf.NUC_DATA_PATH
    with tb.openFile(nuc_data, 'r') as f:
        N = f.root.neutron.eaf_xs.eaf_xs.coldescrs['xsec'].shape[0]
        cond = "(nuc_zz == {0}) & (rxnum == '{1}')".format(nuc, rx)
        rows = [np.array(row['xsec']) for row in \
                f.root.neutron.eaf_xs.eaf_xs.where(cond)]

    if len(rows) == 0:
        # No absportion, return zero-array
        print 'not found'
        sigma_n = np.zeros(N, dtype=float)
    else:
        rows = np.array(rows)
        sigma_n = rows.sum(axis=0)

    return sigma_n


# Make a singleton of the cross-section cache
xs_cache = XSCache()
