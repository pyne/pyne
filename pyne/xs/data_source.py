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
    '''munge the reaction rate name
    '''
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
    
    def __init__(self, src_phi_g=None, dst_group_struct=None, **kwargs):
        """Cross section data source."""
        self._exists = None
        if not self.exists:
            return 
        self._load_group_structure()
        self.dst_group_struct = dst_group_struct
        self.src_phi_g = np.ones(self._src_ngroups, dtype='f8') if src_phi_g is None \
                            else np.asarray(src_phi_g)
        

    @property
    def src_group_struct(self):
        return self._src_group_struct

    @src_group_struct.setter
    def src_group_struct(self, src_group_struct):
        self._src_group_struct = src_group_struct
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

    def discretize(self, nuc, rx, src_phi_g=None, dst_phi_g=None):
        src_phi_g = self.src_phi_g if src_phi_g is None else np.asarray(src_phi_g) 
        src_sigma = self.reaction(nuc, rx)
        dst_sigma = None if src_sigma is None else group_collapse(src_sigma, 
                                                        src_phi_g, dst_phi_g, 
                                                        self._src_to_dst_matrix)
        return dst_sigma

    # Mix-in methods to implement
    def _load_group_structure(self):
        raise NotImplementedError

    @property
    def exists(self):
        raise NotImplementedError

    def reaction(self, nuc, rx):
        raise NotImplementedError


class NullDataSource(DataSource):
    
    def __init__(self, **kwargs):
        """Cross section data source that always returns zeros."""
        super(NullDataSource, self).__init__(**kwargs)

    def _load_group_structure(self):
        """Loads a meaningless bounds array.
        """
        self.src_group_struct = np.array([0.0])

    @property
    def exists(self):
        if self._exists is None:
            self._exists = True
        return self._exists

    def reaction(self, nuc, rx):
        return np.zeros(self.src_ngroups, dtype='f8')

    def discretize(self, nuc, rx, src_phi_g=None, dst_phi_g=None):
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


class CinderDataSource(DataSource):
    
    def __init__(self, **kwargs):
        """Cinder cross section data source."""
        super(CinderDataSource, self).__init__(**kwargs)

    def _load_group_structure(self):
        """Loads and returns the cinder energy bounds array, E_g, 
        from the nuc_data library.
        """
        with tb.openFile(nuc_data, 'r') as f:
            E_g = np.array(f.root.neutron.cinder_xs.E_g)
        self.src_group_struct = E_g

    @property
    def exists(self):
        if self._exists is None:
            with tb.openFile(nuc_data, 'r') as f:
                self._exists = ('/neutron/cinder_xs' in f)
        return self._exists

    def reaction(self, nuc, rx):
        nuc = nucname.zzaaam(nuc)
        rx = _munge_rx(rx)
        f = tb.openFile(nuc_data, 'r')
        node = f.root.neutron.cinder_xs.fission if rx == 'f' else \
               f.root.neutron.cinder_xs.absorption

        # Set query condition
        if rx == 'f':
            cond = 'nuc == {0}'.format(nuc)
        elif rx == 'a':
            cond = "(from_nuc == {0}) & (reaction_type != 'c')".format(nuc)
        elif rx in RX_TYPES:
            cond = "(from_nuc == {0}) & (reaction_type == '{1}')".format(nuc, rx)
        else:
            return None

        rows = [np.array(row['xs']) for row in node.where(cond)]
        f.close()
        if len(rows) == 1:
            rxdata = rows[0]
        elif 1 < len(rows) and rx == 'a':
            rows = np.array(rows)
            rxdata = rows.sum(axis=0)
        else:
            rxdata = None
        return rxdata
