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
    
    def __init__(self, dst_E_g=None, **kwargs):
        """Cross section data source."""
        if not self.exists():
            return 
        self._load_group_structure()
        if dst_E_g is None:
            self.dst_E_g = None
            self.dst_G = 0
            self.src_to_dst_matrix = None
        else:
            self.dst_E_g = np.asarray(dst_E_g)
            self.dst_G = len(dst_E_g) - 1
            self.src_to_dst_matrix = partial_energy_matrix(dst_E_g, self.src_E_g)

    def _load_group_structure(self):
        pass

    def exists(self):
        pass

    def reaction(self, nuc, rx):
        pass

    def collapse(self, nuc, rx, src_phi_g=None, dst_phi_g=None):
        src_phi_g = np.asarray(src_phi_g) if src_phi_g is not None else \
                    np.ones(self.src_G, dtype='float64')
        src_sigma = self.reaction(nuc, rx)
        dst_sigma = group_collapse(src_sigma, src_phi_g, dst_phi_g, 
                                   self.src_to_dst_matrix)
        return dst_sigma


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
        self.src_G = len(E_g) - 1
        self.src_E_g = E_g

    def exists(self):
        with tb.openFile(nuc_data, 'r') as f:
            rtn = ('/neutron/cinder_xs' in f)
        return rtn

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
