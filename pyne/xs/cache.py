"""This module provides a cross section cache which automatically extracts 
cross-sections from provided nuclear data sets."""
from itertools import product
from collections import MutableMapping

import numpy as np
import tables as tb

from pyne import nucname
from pyne.pyne_config import pyne_conf
from pyne.xs.models import partial_energy_matrix, phi_g
from pyne.xs import data_source


def _same_arr_or_none(a, b): 
    if a is None or b is None:
        return a is b
    else:
        return (len(a) == len(b)) and (a == b).all()

def _valid_group_struct(E_g):
    if E_g is None:
        return None
    E_g = np.asarray(E_g, dtype='f8')
    if E_g[0] < E_g[-1]:
        E_g = E_g[::-1]
    return E_g


###############################################################################
### Set up a cross-section cache so the same data isn't loaded repetitively ###
###############################################################################

class XSCache(MutableMapping):
    """A lightweight multigroup cross section cache based off of python 
    dictionaries. This relies on a list of cross section data source from which
    discretized group constants may be computed from raw, underlying data.  
    Normally this requires that some cross section data be built into nuc_data.
    A default instance of this class is provided (pyne.xs.cache.xs_cache).

    Parameters
    ----------
    group_struct : array-like of floats, optional
        Energy group structure E_g [MeV] from highest-to-lowest energy, length G+1.
        If the group structure is not present in the cache or is None, all cross 
        sections pulled from the sources will not be discretized or collapsed.
    data_source_classes : list of DataSource classes, optional
        Sequence of DataSource classes (not instances!) from which to grab cross 
        section data. Data from a source earlier in the sequence (eg, index 1)
        will take precednce over data later in the sequence (eg, index 5).

    """

    def __init__(self, group_struct=None, 
                 data_source_classes=(data_source.CinderDataSource,
                                      data_source.SimpleDataSource,
                                      data_source.NullDataSource,)):
        self._cache = {}
        self.data_sources = []
        for cls in data_source_classes:
            ds = cls(dst_group_struct=group_struct)
            if ds.exists:
                self.data_sources.append(ds)
        self._cache['E_g'] = _valid_group_struct(group_struct) 
        self._cache['phi_g'] = None

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
        if (key not in self._cache) and not isinstance(key, basestring):
            E_g = self._cache['E_g']
            if E_g is None:
                for ds in self.data_sources:
                    xsdata = ds.reaction(*key)
                    if xsdata is not None:
                        self._cache[key] = xsdata
                        break
            else:
                kw = dict(zip(['nuc', 'rx', 'temp'], key))
                kw['dst_phi_g'] = self._cache['phi_g']
                for ds in self.data_sources:
                    xsdata = ds.discretize(**kw)
                    if xsdata is not None:
                        self._cache[key] = xsdata
                        break            
        # Return the value requested
        return self._cache[key]


    def __setitem__(self, key, value):
        """Key setting via custom cache functionality."""
        # Set the E_g
        if (key == 'E_g'):
            value = _valid_group_struct(value)
            cache_value = self._cache['E_g']
            if _same_arr_or_none(value, cache_value):
                return 
            self.clear()
            self._cache['phi_g'] = None
            for ds in self.data_sources:
                ds.dst_group_struct = value
        elif (key == 'phi_g'):
            value = value if value is None else np.asarray(value, dtype='f8')
            cache_value = self._cache['phi_g']
            if _same_arr_or_none(value, cache_value):
                return
            E_g = self._cache['E_g']
            if len(value) + 1 == len(E_g):
                self.clear()
            else:
                raise ValueError("phi_g does not match existing group structure E_g!")
        # Set the value normally
        self._cache[key] = value

    def clear(self):
        """Clears the cache, retaining E_g and phi_g."""
        E_g, phi_g = self._cache['E_g'], self._cache['phi_g'] 
        self._cache.clear()
        self._cache['E_g'], self._cache['phi_g'] = E_g, phi_g
        


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

EAF_RX_MAP = { #TODO
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
    """Grabs group-wise cross section values from a given nuclide and reaction from the EAF format activation library.

    Parameters
    ----------
    nuc : int
        Nuclide in zzaaam form.
    rx : str 
        Reaction MT # in nnnm form.
        OR: (eventually)
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


