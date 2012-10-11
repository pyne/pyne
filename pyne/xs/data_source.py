"""Cross section library data source interfaces.
"""
import numpy as np
import tables as tb

from pyne import nuc_data
from pyne import nucname
from pyne.xs.models import partial_energy_matrix, group_collapse

RX_TYPES = set(["", 'np *', 'a  *', 'h  *', '2p *', '3n *', 'd  *', 'np/d',
                'na', '*', 'nd', 'g  *', '3n', 'np', 'nt', 't', 'nt *',
                'x  *', '4n *', 'na *', 'nd *', 't  *', 'a', 'c', '2p', 'd',
                'g', 'h', 'n', '4n', 'p', 'n  *', '2a', '2n *', 'x', '2n',
                'nh *', 'p  *', 'f'])

RX_TYPES_MAP = {
    'neutron': 'n',
    'gamma': 'g',
    'alpha': 'a',
    'proton': 'p',
    'trit': 't',
    'triton': 't',
    'deut': 'd',
    'deuteron': 'd',
    }

def _munge_rx(rx):
    """Munge the reaction rate name."""
    if rx not in RX_TYPES:
        while any([(key in rx) for key in RX_TYPES_MAP]):
            for key, value in RX_TYPES_MAP.items():
                rx = rx.replace(key, value)

    if '_x' in rx:
        if len(rx) == 3:
            rx = rx.replace('_x', '  *')
        else:
            rx = rx.replace('_x', ' *')

    if rx not in RX_TYPES:
        msg = "the reaction '{rx}' is not valid.".format(rx=rx)
        raise IndexError(msg)
    return rx


class DataSource(object):
    """Base cross section data source.

    This is an abstract class which provides default functionality when subclassed
    and certain methods are overridden.  A data source should know how to find cross 
    section information for a given nuclide, reaction type, and temperature, if 
    available.  If such data is not present for this source, then the source should
    return None.

    Furthermore, a data source must be able to declare whether the data is present
    on the users system.

    Finally, data sources distinguish between group structure information coming
    from or for the source itself (src) and optional user-defined destination (dst) 
    group structure information.  This allows the data source to define how it wishes
    to discretize its data to the custom destination form.

    The following methods must be overridden in all DataSource subclasses:

    .. code-block:: python

        @property
        def exists(self):
            # Is this DataSource available on the user's system?
            return (True or False)

        def _load_group_structure(self):
            # Sets the source group structure, E_g, native to this DataSource
            ...
            self.src_group_struct = E_g (array-like)

        def _load_reaction(self, nuc, rx, temp=300.0):
            # Returns the rection channel for a nuclide at a given temperature
            # or returns None, if this unavailable in this DataSource.
            ...
            return rxdata (ndarray of floats, length self.src_ngroups, or None)

    Note that non-multigroup data sources should also override the discretize()
    method.  Other methods and properties may also need to be overriden depending
    on the data source at hand.

    All data sources may be used independenly or in conjunction with a cross
    section cache instance.

    Parameters
    ----------
    src_phi_g : array-like, optional
        Group fluxes which must match the group structure for this data source.
    dst_group_struct : array-like, optional
        The energy group structure [MeV] of the destination cross sections.  Used when 
        discretizing cross sections from this source.

    """
    
    def __init__(self, src_phi_g=None, dst_group_struct=None, **kwargs):
        self._exists = None
        if not self.exists:
            return 
        self.rxcache = {}
        self._load_group_structure()
        self.dst_group_struct = dst_group_struct
        self.src_phi_g = np.ones(self._src_ngroups, dtype='f8') if src_phi_g is None \
                            else np.asarray(src_phi_g)
        

    @property
    def src_group_struct(self):
        return self._src_group_struct

    @src_group_struct.setter
    def src_group_struct(self, src_group_struct):
        self._src_group_struct = np.asarray(src_group_struct, dtype='f8')
        self._src_ngroups = len(src_group_struct) - 1        

    @property
    def src_ngroups(self):
        return self._src_ngroups

    @property
    def dst_group_struct(self):
        return self._dst_group_struct

    @dst_group_struct.setter
    def dst_group_struct(self, dst_group_struct):
        if dst_group_struct is None:
            self._dst_group_struct = None
            self._dst_ngroups = 0
            self._src_to_dst_matrix = None
        else:
            self._dst_group_struct = np.asarray(dst_group_struct)
            self._dst_ngroups = len(dst_group_struct) - 1
            self._src_to_dst_matrix = partial_energy_matrix(dst_group_struct, 
                                                            self._src_group_struct)

    @property
    def dst_ngroups(self):
        return self._dst_ngroups

    @property
    def src_to_dst_matrix(self):
        return self._src_to_dst_matrix

    def reaction(self, nuc, rx, temp=300.0):
        """Gets the cross section data for this reaction channel either directly from
        the data source or from the rxcache.

        Parameters
        ----------
        nuc : int or str
            A nuclide.
        rx : int or str
            Reaction key ('gamma', 'alpha', 'p', etc.) or MT number.
        temp : float, optional
            Temperature [K] of material, defaults to 300.0.

        Returns
        -------
        rxdata : ndarry
            Source cross section data, length src_ngroups.

        """
        rxkey = (nuc, rx, temp)
        if rxkey not in self.rxcache:
            self.rxcache[rxkey] = self._load_reaction(nuc, rx, temp)
        return self.rxcache[rxkey]

    def discretize(self, nuc, rx, temp=300.0, src_phi_g=None, dst_phi_g=None):
        """Discretizes the reaction channel from the source group structure to that 
        of the destination weighted by the group fluxes.  This implemenation is only
        valid for multi-group data sources.  Non-multigroup data source should also
        override this method.

        Parameters
        ----------
        nuc : int or str
            A nuclide.
        rx : int or str
            Reaction key ('gamma', 'alpha', 'p', etc.) or MT number.
        temp : float, optional
            Temperature [K] of material, defaults to 300.0.
        src_phi_g : array-like, optional
            Group fluxes for this data source, length src_ngroups.
        dst_phi_g : array-like, optional
            Group fluxes for the destiniation structure, length dst_ngroups.

        Returns
        -------
        dst_sigma : ndarry
            Destination cross section data, length dst_ngroups.

        """
        src_phi_g = self.src_phi_g if src_phi_g is None else np.asarray(src_phi_g) 
        src_sigma = self.reaction(nuc, rx, temp)
        dst_sigma = None if src_sigma is None else group_collapse(src_sigma, 
                                                        src_phi_g, dst_phi_g, 
                                                        self._src_to_dst_matrix)
        return dst_sigma

    # Mix-in methods to implement
    @property
    def exists(self):
        raise NotImplementedError

    def _load_group_structure(self):
        raise NotImplementedError

    def _load_reaction(self, nuc, rx, temp=300.0):
        raise NotImplementedError


class NullDataSource(DataSource):
    """Cross section data source that always exists and always returns zeros.

    Parameters
    ----------
    kwargs : optional
        Keyword arguments to be sent to base class.

    """
    
    def __init__(self, **kwargs):
        super(NullDataSource, self).__init__(**kwargs)

    def _load_group_structure(self):
        """Loads a meaningless bounds array."""
        self.src_group_struct = np.array([0.0])

    @property
    def exists(self):
        if self._exists is None:
            self._exists = True
        return self._exists

    def _load_reaction(self, nuc, rx, temp=300.0):
        return np.zeros(self.src_ngroups, dtype='f8')

    def discretize(self, nuc, rx, temp=300.0, src_phi_g=None, dst_phi_g=None):
        """Returns zeros."""
        return np.zeros(self.dst_ngroups, dtype='f8')

    @property
    def dst_group_struct(self):
        return self._dst_group_struct

    @dst_group_struct.setter
    def dst_group_struct(self, dst_group_struct):
        if dst_group_struct is None:
            self._dst_group_struct = None
            self._dst_ngroups = 0
        else:
            self._dst_group_struct = np.asarray(dst_group_struct)
            self._dst_ngroups = len(dst_group_struct) - 1
        self._src_to_dst_matrix = None


class SimpleDataSource(DataSource):
    """Simple cross section data source based off of KAERI data.

    Parameters
    ----------
    kwargs : optional
        Keyword arguments to be sent to base class.

    """
    _rx_avail = set(['t', 's', 'e', 'i', 'a', 'gamma', 'f', 'alpha', 'proton', 
                     'deut', 'trit', '2n', '3n', '4n'])

    def __init__(self, **kwargs):
        super(SimpleDataSource, self).__init__(**kwargs)

    @property
    def exists(self):
        if self._exists is None:
            with tb.openFile(nuc_data, 'r') as f:
                self._exists = ('/neutron/simple_xs' in f)
        return self._exists

    def _load_group_structure(self):
        """Sets the simple energy bounds array, E_g."""
        self.src_group_struct = np.array([14.0, 1.0, 2.53E-8, 0.0], dtype='float64')

    def _load_reaction(self, nuc, rx, temp=300.0):
        nuc = nucname.zzaaam(nuc)
        if rx not in self._rx_avail:
            return None
        cond = "nuc == {0}".format(nuc)
        sig = 'sigma_' + rx
        with tb.openFile(nuc_data, 'r') as f:
            simple_xs = f.root.neutron.simple_xs
            fteen = [row[sig] for row in simple_xs.fourteen_MeV.where(cond)]
            fissn = [row[sig] for row in simple_xs.fission_spectrum_ave.where(cond)]
            therm = [row[sig] for row in simple_xs.thermal.where(cond)]
            if 0 == len(therm):
                therm = [row[sig] for row in simple_xs.thermal_maxwell_ave.where(cond)]
            if 0 == len(fteen) and 0 == len(fissn) and 0 == len(therm):
                rxdata = None
            else:
                rxdata = np.array([fteen[0], fissn[0], therm[0]], dtype='float64')
        return rxdata

    def discretize(self, nuc, rx, temp=300.0, src_phi_g=None, dst_phi_g=None):
        """Discretizes the reaction channel from simple group structure to that 
        of the destination weighted by the group fluxes.  Since the simple data 
        source consists of only thermal (2.53E-8 MeV), fission (1 MeV), and 14 MeV
        data points, the following piecewise functional form is assumed:

        .. math::  

            \sigma(E) = \sigma(2.53E-8) \sqrt{\frac{2.53E-8}{E}} 
            \sigma(E) = \frac{\sigma(14) - \sigma(1)}{14 - 1} (E - 1) + \sigma(1) 

        Parameters
        ----------
        nuc : int or str
            A nuclide.
        rx : int or str
            Reaction key ('gamma', 'alpha', 'p', etc.) or MT number.
        temp : float, optional
            Temperature [K] of material, defaults to 300.0.
        src_phi_g : array-like, optional
            IGNORED!!!  Included for API compatability
        dst_phi_g : array-like, optional
            Group fluxes for the destiniation structure, length dst_ngroups.

        Returns
        -------
        dst_sigma : ndarry
            Destination cross section data, length dst_ngroups.

        """
        src_phi_g = self.src_phi_g if src_phi_g is None else np.asarray(src_phi_g) 
        src_sigma = self.reaction(nuc, rx, temp)
        if src_sigma is None:
            return None
        # not the most efficient, but data sizes should be smallish
        center_g = self._dst_centers
        dst_sigma = (src_sigma[2] * np.sqrt(2.53E-8)) / np.sqrt(center_g)
        dst_fissn = ((src_sigma[0] - src_sigma[1])/13.0) * (center_g - 1.0) + \
                                                                        src_sigma[1]
        mask = (dst_sigma < dst_fissn)
        dst_sigma[mask] = dst_fissn[mask]
        if dst_phi_g is not None:
            dst_sigma = (dst_sigma * dst_phi_g) / dst_phi_g.sum()
        return dst_sigma

    @property
    def dst_group_struct(self):
        return self._dst_group_struct

    @dst_group_struct.setter
    def dst_group_struct(self, dst_group_struct):
        if dst_group_struct is None:
            self._dst_group_struct = None
            self._dst_centers = None
            self._dst_ngroups = 0
        else:
            self._dst_group_struct = np.asarray(dst_group_struct)
            self._dst_centers = (self._dst_group_struct[1:] + 
                                 self._dst_group_struct[:-1])/2.0
            self._dst_ngroups = len(dst_group_struct) - 1
        self._src_to_dst_matrix = None


class CinderDataSource(DataSource):
    """Cinder cross section data source. The relevant cinder cross section data must
    be present in the nuc_data for this data source to exist.

    Parameters
    ----------
    kwargs : optional
        Keyword arguments to be sent to base class.

    """
    
    def __init__(self, **kwargs):
        super(CinderDataSource, self).__init__(**kwargs)

    def _load_group_structure(self):
        """Loads the cinder energy bounds array, E_g, from nuc_data."""
        with tb.openFile(nuc_data, 'r') as f:
            E_g = np.array(f.root.neutron.cinder_xs.E_g)
        self.src_group_struct = E_g

    @property
    def exists(self):
        if self._exists is None:
            with tb.openFile(nuc_data, 'r') as f:
                self._exists = ('/neutron/cinder_xs' in f)
        return self._exists

    def _load_reaction(self, nuc, rx, temp=300.0):
        nuc = nucname.zzaaam(nuc)
        rx = _munge_rx(rx)

        # Set query condition
        if rx == 'f':
            cond = 'nuc == {0}'.format(nuc)
        elif rx in RX_TYPES:
            cond = "(from_nuc == {0}) & (reaction_type == '{1}')".format(nuc, rx)
        else:
            return None

        with tb.openFile(nuc_data, 'r') as f:
            node = f.root.neutron.cinder_xs.fission if rx == 'f' else \
                   f.root.neutron.cinder_xs.absorption
            rows = [np.array(row['xs']) for row in node.where(cond)]

        if 1 == len(rows):
            rxdata = rows[0]
        elif 1 < len(rows):
            rows = np.array(rows)
            rxdata = rows.sum(axis=0)
        elif 0 == len(rows) and (rx == 'a'):
            # in case absorption doesn't exist, we compute it
            fdata = self._load_reaction(nuc, 'f')
            cond = "(from_nuc == {0}) & (reaction_type != 'c')".format(nuc)
            with tb.openFile(nuc_data, 'r') as f:
                node = f.root.neutron.cinder_xs.absorption
                rows = np.array([row['xs'] for row in node.where(cond)])
            if 0 == len(rows) and fdata is None:
                rxdata = None
            else:
                rxdata = rows.sum(axis=0)
                if fdata is not None:
                    rxdata += fdata
        else:
            rxdata = None
        return rxdata



class EAFDataSource(DataSource):
    """European Activation File cross section data source.  The relevant EAF
    cross section data must be present in the nuc-data for this data source to exist.

    Parameters
    ----------

    """

    def __init__(self):
        super(EAFDataSource, self).__init__()

    def _load_group_structure(self):
        """ """
        with tb.openFile(nuc_data, 'r') as f:
            E_g = np.array(f.root.neutron.eaf_xs.E_g
