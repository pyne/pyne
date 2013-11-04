"""Cross section library data source interfaces.
"""
import os
import StringIO

import numpy as np
import tables as tb

from pyne import nuc_data
from pyne import nucname
from pyne import rxname
from pyne.xs.models import partial_energy_matrix, group_collapse
from pyne.endf import Library

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

    The following methods may be overridden in DataSource subclasses as a potential
    optimization:

    .. code-block:: python

        def load(self, temp=300.0):
            # loads the entire data source into memory.  This prevents 
            # excessive queries to disk.  This does not return anything.
            pass

        # a class attribute that specificies whether the data source uses 
        # temperature information.  If the data source does not use such 
        # info, the keys in the rxcache dictionary are shortend to just the 
        # (nuc, rx) pair.  Accessing (nuc, rx, temp) data will defer to the 
        # (nuc, rx) pair.
        _USES_TEMP = True

    Note that non-multigroup data sources should also override the discretize()
    method.  Other methods and properties may also need to be overriden depending
    on the data source at hand.

    All data sources may be used independently or in conjunction with a cross
    section cache instance.

    Parameters
    ----------
    src_phi_g : array-like, optional
        Group fluxes which must match the group structure for this data source.
    dst_group_struct : array-like, optional
        The energy group structure [MeV] of the destination cross sections.
        Used when discretizing cross sections from this source.

    """

    def __init__(self, src_phi_g=None, dst_group_struct=None, **kwargs):
        self._exists = None
        if not self.exists:
            return
        self.rxcache = {}
        self.fullyloaded = False
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
        """Gets the cross section data for this reaction channel either directly
        from the data source or from the rxcache.

        Parameters
        ----------
        nuc : int or str
            A nuclide.
        rx : int or str
            Reaction id or name.
        temp : float, optional
            Temperature [K] of material, defaults to 300.0.

        Returns
        -------
        rxdata : ndarray
            Source cross section data, length src_ngroups.

        """
        nuc = nucname.id(nuc)
        rx = rxname.id(rx)
        rxkey = (nuc, rx, temp) if self._USES_TEMP else (nuc, rx)
        if rxkey not in self.rxcache:
            self.rxcache[rxkey] = None if self.fullyloaded \
                                       else self._load_reaction(nuc, rx, temp)
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
            Reaction id or name.
        temp : float, optional
            Temperature [K] of material, defaults to 300.0.
        src_phi_g : array-like, optional
            Group fluxes for this data source, length src_ngroups.
        dst_phi_g : array-like, optional
            Group fluxes for the destiniation structure, length dst_ngroups.

        Returns
        -------
        dst_sigma : ndarray
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

    # Optional mix-in methods to implement
    def load(self, temp=300.0):
        pass

    _USES_TEMP = True

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

    _USES_TEMP = False

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
    """Simple cross section data source based off of KAERI data.  This data source
    does not use material temperature information.

    Parameters
    ----------
    kwargs : optional
        Keyword arguments to be sent to base class.

    """
    _rx_avail = {rxname.id('total'): 't', 
                 rxname.id('scattering'): 's', 
                 rxname.id('elastic'): 'e', 
                 rxname.id('inelastic'): 'i', 
                 rxname.id('absorption'): 'a', 
                 rxname.id('gamma'): 'gamma', 
                 rxname.id('fission'): 'f', 
                 rxname.id('alpha'): 'alpha', 
                 rxname.id('proton'): 'proton', 
                 rxname.id('deut'): 'deut', 
                 rxname.id('trit'): 'trit', 
                 rxname.id('z_2n'): '2n', 
                 rxname.id('z_3n'): '3n', 
                 rxname.id('z_4n'): '4n'}

    def __init__(self, **kwargs):
        super(SimpleDataSource, self).__init__(**kwargs)

    @property
    def exists(self):
        if self._exists is None:
            with tb.openFile(nuc_data, 'r') as f:
                self._exists = ('/neutron/simple_xs' in f)
        return self._exists

    _USES_TEMP = False

    def _load_group_structure(self):
        """Sets the simple energy bounds array, E_g."""
        self.src_group_struct = np.array([14.0, 1.0, 2.53E-8, 0.0], dtype='float64')

    def _load_reaction(self, nuc, rx, temp=300.0):
        if rx not in self._rx_avail:
            return None
        cond = "nuc == {0}".format(nuc)
        sig = 'sigma_' + self._rx_avail[rx]
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

            \\sigma(E) = \\sigma(2.53E-8) \\sqrt{\\frac{2.53E-8}{E}}
            \\sigma(E) = \\frac{\sigma(14) - \\sigma(1)}{14 - 1} (E - 1) + \\sigma(1)

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
        dst_sigma : ndarray
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
    be present in the nuc_data for this data source to exist.  This data source does
    not use material temperature information.

    Parameters
    ----------
    kwargs : optional
        Keyword arguments to be sent to base class.

    """
    # 'h' stands for helion or 'He3'
    _rx_avail = {rxname.id('np_1'): 'np *',
                 rxname.id('a_1'): 'a  *', 
                 rxname.id('He3_1'): 'h  *', 
                 rxname.id('z_2p_1'): '2p *', 
                 rxname.id('z_3n_1'): '3n *', 
                 rxname.id('d_1'): 'd  *',
                 rxname.id('npd'): 'np/d',
                 rxname.id('na'): 'na', 
                 rxname.id('excited'): '*', 
                 rxname.id('nd'): 'nd', 
                 rxname.id('gamma_1'): 'g  *', 
                 rxname.id('z_3n'): '3n', 
                 rxname.id('np'): 'np', 
                 rxname.id('nt'): 'nt', 
                 rxname.id('t'): 't', 
                 rxname.id('nt_1'): 'nt *',
                 rxname.id('z_4n_1'): '4n *', 
                 rxname.id('na_1'): 'na *', 
                 rxname.id('nd_1'): 'nd *', 
                 rxname.id('t_1'): 't  *', 
                 rxname.id('a'): 'a', 
                 rxname.id('z_2p'): '2p', 
                 rxname.id('d'): 'd',
                 rxname.id('gamma'): 'g', 
                 rxname.id('He3'): 'h', 
                 rxname.id('n'): 'n', 
                 rxname.id('z_4n'): '4n', 
                 rxname.id('p'): 'p', 
                 rxname.id('n_1'): 'n  *', 
                 rxname.id('z_2a'): '2a', 
                 rxname.id('z_2n_1'): '2n *', 
                 rxname.id('z_2n'): '2n',
                 rxname.id('nHe3_1'): 'nh *', 
                 rxname.id('p_1'): 'p  *', 
                 # not real or unique absorption reactions
                 #rxname.id(''): "", 
                 #rxname.id(''): 'x', 
                 #rxname.id(''): 'x  *',
                 #rxname.id(''): 'c', 
                 #rxname.id('fission'): 'f',
                 }

    _USES_TEMP = False

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
        fissrx = rxname.id('fission')
        absrx = rxname.id('absorption')

        # Set query condition
        if rx in self._rx_avail:
            cond = "(from_nuc == {0}) & (reaction_type == '{1}')"
            cond = cond.format(nuc, self._rx_avail[rx])
        elif rx == fissrx:
            cond = 'nuc == {0}'.format(nuc)
        elif rx == absrx:
            cond = "(from_nuc == {0}) & (reaction_type != 'c')".format(nuc)
        else:
            return None

        # read & collapse data
        with tb.openFile(nuc_data, 'r') as f:
            node = f.root.neutron.cinder_xs.fission if rx == fissrx else \
                   f.root.neutron.cinder_xs.absorption
            rows = [np.array(row['xs']) for row in node.where(cond)]

        if 1 == len(rows):
            rxdata = rows[0]
        elif 1 < len(rows):
            rows = np.array(rows)
            rxdata = rows.sum(axis=0)
        else:
            rxdata = None

        # add fission data to absorption
        if rx == absrx:
            fdata = self._load_reaction(nuc, fissrx)
            if fdata is not None:
                rxdata = fdata if rxdata is None else rxdata + fdata
        return rxdata

class EAFDataSource(DataSource):
    """European Activation File cross section data source.  The relevant EAF
    cross section data must be present in the nuc-data for this data source
    to exist.

    Parameters
    ----------
    kwargs : optional
        Keyword arguments to be sent to base class.

    Notes
    -----
    EAF data does not use temperature information.
    """

    # MT#s included in the EAF data
    _rx_avail = {rxname.id('gamma'): '1020', 
                 rxname.id('gamma_1'): '1021', 
                 rxname.id('gamma_2'): '1022', 
                 rxname.id('p'): '1030', 
                 rxname.id('p_1'): '1031',
                 rxname.id('p_2'): '1032',
                 rxname.id('d'): '1040', 
                 rxname.id('d_1'): '1041', 
                 rxname.id('d_2'): '1042',
                 rxname.id('t'): '1050',
                 rxname.id('t_1'): '1051',
                 rxname.id('t_2'): '1052',
                 rxname.id('He3'): '1060',
                 rxname.id('He3_1'): '1061',
                 rxname.id('He3_2'): '1062',
                 rxname.id('a'): '1070',
                 rxname.id('a_1'): '1071',
                 rxname.id('a_2'): '1072',
                 rxname.id('z_2a'): '1080',
                 rxname.id('z_2p'): '1110',
                 rxname.id('z_2p_1'): '1111',
                 rxname.id('z_2p_2'): '1112',
                 rxname.id('z_2n'): '160',
                 rxname.id('z_2n_1'): '161',
                 rxname.id('z_2n_2'): '162',
                 rxname.id('z_3n'): '170',
                 rxname.id('z_3n_1'): '171',
                 rxname.id('z_3n_2'): '172',
                 rxname.id('fission'): '180',
                 rxname.id('na'): '220',
                 rxname.id('na_1'): '221',
                 rxname.id('na_2'): '222',
                 rxname.id('z_2na'): '240',
                 rxname.id('np'): '280',
                 rxname.id('np_1'): '281',
                 rxname.id('np_2'): '282',
                 rxname.id('n2a'): '290',
                 rxname.id('nd'): '320',
                 rxname.id('nd_1'): '321',
                 rxname.id('nd_2'): '322',
                 rxname.id('nt'): '330',
                 rxname.id('nt_1'): '331',
                 rxname.id('nt_2'): '332',
                 rxname.id('nHe3'): '340',
                 rxname.id('nHe3_1'): '341',
                 rxname.id('nHe3_2'): '342',
                 rxname.id('z_4n'): '370',
                 rxname.id('z_4n_1'): '371',
                 rxname.id('n'): '40',
                 rxname.id('n_1'): '41',
                 rxname.id('n_2'): '42',
                 rxname.id('z_3np'): '420',
                 }

    _avail_rx = dict([_[::-1] for _ in _rx_avail.items()])

    _USES_TEMP = False

    def __init__(self, **kwargs):
        super(EAFDataSource, self).__init__(**kwargs)

    def _load_group_structure(self):
        """Loads the EAF energy bounds array, E_g, from nuc_data."""
        with tb.openFile(nuc_data, 'r') as f:
            E_g = np.array(f.root.neutron.eaf_xs.E_g)
        self.src_group_struct = E_g

    @property
    def exists(self):
        if self._exists is None:
            with tb.openFile(nuc_data, 'r') as f:
                self._exists = ('/neutron/eaf_xs' in f)
        return self._exists

    def _load_reaction(self, nuc, rx, temp=300.0):
        """Loads reaction specific data for EAF.

        Parameters
        ----------
        nuc : int
            Nuclide id.
        rx : int 
            Reaction id.
        temp : float, optional
            The material temperature
        
        Notes
        -----
        EAF data does not use temperature information (temp).

        """
        absrx = rxname.id('absorption')

        if rx in self._rx_avail:
            cond = "(nuc_zz == {0}) & (rxnum == '{1}')".format(nuc, self._rx_avail[rx])
        elif rx == absrx:
            cond = "(nuc_zz == {0})".format(nuc)            
        else:
            return None

        # Grab data
        with tb.openFile(nuc_data, 'r') as f:
            node = f.root.neutron.eaf_xs.eaf_xs
            rows = node.readWhere(cond)
            #rows = [np.array(row['xs']) for row in node.where(cond)]

        if len(rows) == 0:
            rxdata = None
        elif 1 < len(rows):
            xss = rows['xs']
            rxnums = rows['rxnum']
            for rxnum, xs in zip(rxnums, xss):
                self.rxcache[nuc, self._avail_rx[rxnum]] = xs
            rxdata = xss.sum(axis=0)
        else:
            rxdata = rows[0]['xs']

        return rxdata

    def load(self, temp=300.0):
        """Loads all EAF into memory.

        Parameters
        ----------
        temp : float, optional
            The material temperature
        
        Notes
        -----
        EAF data does not use temperature information (temp).

        """
        rxcache = self.rxcache
        avail_rx = self._avail_rx
        absrx = rxname.id('absorption')
        with tb.openFile(nuc_data, 'r') as f:
            node = f.root.neutron.eaf_xs.eaf_xs
            for row in node:
                nuc = row['nuc_zz']
                rx = avail_rx[row['rxnum']]
                xs = row['xs']
                rxcache[nuc, rx] = xs
                abskey = (nuc, absrx)
                rxcache[abskey] = xs + rxcache.get(abskey, 0.0)            
        self.fullyloaded = True


class ENDFDataSource(DataSource):
    """Evaluated Nuclear Data File cross section data source.  The ENDF file
    must exist for this data source to exist.

    Parameters
    ----------
    f : string, file handle
        Path to ENDF file, or ENDF file itself.
    kwargs : optional
        Keyword arguments to be sent to base class.

    Notes
    -----
    """

    def __init__(self, fh, src_phi_g=None, dst_group_struct=None, **kwargs):
        self.fh = fh
        self._exists = None
        if not self.exists:
            raise ValueError
        else:
            self.library = Library(fh)
        self.rxcache = {}
        self.dst_group_struct = dst_group_struct
        self._src_phi_g = src_phi_g

    @property
    def exists(self):
        if self._exists is None:
            if isinstance(self.fh, basestring):
                self._exists = os.path.isfile(fh)
            else:
                self._exists = (isinstance(self.fh, file) or \
                                isinstance(self.fh, StringIO.StringIO))
        return self._exists

    def _load_group_structure(self, nuc, rx, nuc_i=None):
        """Loads the group structure from ENDF file."""
        self.library._read_res(nuc)
        mt = rxname.mt(rx)
        rx_data = self.rxcache[nuc, rx, nuc_i]
        xsdata = self.library.get_xs(nuc, rx, nuc_i)[0]
        intpoints = xsdata['intpoints']#[::-1]
        Eint = xsdata['Eint']
        E_g = []
        for i in range(len(intpoints)):
            if not i:
                low_Eint = 0
            else:
                low_Eint = intpoints[i-1]
            high_Eint = intpoints[i]
            E_g.append(Eint[low_Eint:high_Eint])
        rx_data['src_group_struct'] = E_g
        rx_data['src_phi_g'] = np.ones(len(E_g), dtype='f8') \
            if rx_data['_src_phi_g'] is None \
            else np.asarray(rx_data['src_phi_g'])

    def _load_reaction(self, nuc, rx, nuc_i, src_phi_g=None, temp=300.0):
        """Note: EAF data does not use temperature information (temp)

        Parameters
        ----------
        nuc : int
            Nuclide in id form.
        rx : int or str
            Reaction MT # in nnnm form.
            OR:
            Reaction key: 'gamma', 'alpha', 'p', etc.
        nuc_i: : int
            Isotope in id form (optional). Default is None.

        Returns
        -------
        rxdata: ndarray of floats, len ngroups
        """
        nuc = nucname.id(nuc)
        # Munging the rx to an MT#
        try:
            rx = int(rx)
        except ValueError:
            rx = rxname.mt(rx)
        # Check if usable rx #
        if rx is None:
            return None
        # Grab data
        if (nuc, rx, nuc_i) in self.rxcache:
            rxdict = self.rxcache[nuc, rx, nuc_i]
        else:
            if nuc_i not in self.library.structure[nuc]['data']:
                self.library._read_res(nuc)
            rxdict = self.library.get_xs(nuc, rx, nuc_i)[0]
            rxdict['_src_phi_g'] = src_phi_g
            self.rxcache[nuc, rx, nuc_i] = rxdict
        self._load_group_structure(nuc, rx, nuc_i)
        rxdata = rxdict['Eint']
        return rxdata

    def discretize(self, nuc, rx, nuc_i, temp=300.0, src_phi_g=None,
                   dst_phi_g=None):
        """Discretizes the reaction channel.

        Parameters
        ----------
        nuc : int
            Nuclide to discretize.
        rx : int or str
            Reaction to discretize.
        nuc_i : int
            Isotope to discretize.
        temp : float
            Temperature - not used in this, but preserved for API compatibility
        src_phi_g : array-like
            Source group flux - not used in this, but preserved for API
            compatibility
        dst_phi_g : array-like
            Destination group flux - not used in this, but preserved for API
            compatibility

        Returns
        -------
        dst_sigma : array
            An array with the group cross-sections in order of decreasing energy.
        """
        # Munging the rx to an MT#
        try:
            rx = int(rx)
        except ValueError:
            rx = rxname.mt(rx)
        # Check if usable rx #
        if rx is None:
            return None
        self._load_group_structure(nuc, rx, nuc_i)
        self._load_reaction(nuc, rx, nuc_i)
        rxdata = self.rxcache[nuc, rx, nuc_i]
        intpoints = rxdata['intpoints']
        intschemes = rxdata['intschemes']
        Eints = rxdata['Eint']
        E_g = rxdata['src_group_struct']
        xs = rxdata['xs']
        xs_all = []
        dst_sigma = []
        for i in range(len(intpoints)):
            if not i:
                low_xs = 0
            else:
                low_xs = intpoints[i-1]
            high_xs = intpoints[i]
            xs_all.append(xs[low_xs:high_xs])
        for i in range(len(E_g)):
            intscheme = intschemes[i]
            Eints = E_g[i]
            xs = xs_all[i]
            dst_sigma.append(self.library.integrate_tab_range(intscheme, Eints, xs))
        rxdata['dst_sigma'] = np.asarray(dst_sigma)[::-1]
        return rxdata['dst_sigma']
