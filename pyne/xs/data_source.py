"""Cross section library data source interfaces.
"""
from __future__ import division

import os
import io
import sys
from pyne.utils import QA_warn

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import numpy as np
import tables as tb

try:
    # openmc is not a hard dependency of pyne
    from openmc import statepoint
except ImportError:
    statepoint = None

from pyne import nuc_data
from pyne import nucname
from pyne import openmc_utils
from pyne import rxname
from pyne import endf
from pyne import bins
from pyne import ace
from pyne.data import MeV_per_K
from pyne.xs.models import partial_energy_matrix, group_collapse, same_arr_or_none

QA_warn(__name__)

IO_TYPES = (io.IOBase, StringIO)

py3k = False
if sys.version_info[0] > 2:
    basestring = str
    py3k = True


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
        self.src_phi_g = (
            np.ones(self._src_ngroups, dtype="f8")
            if src_phi_g is None
            else np.asarray(src_phi_g)
        )
        self.atom_dens = {}
        self.slf_shld_wgts = {}

    @property
    def src_group_struct(self):
        return self._src_group_struct

    @src_group_struct.setter
    def src_group_struct(self, src_group_struct):
        self._src_group_struct = np.asarray(src_group_struct, dtype="f8")
        self._src_ngroups = len(src_group_struct) - 1
        self.rxcache.clear()

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
        elif same_arr_or_none(dst_group_struct, self._dst_group_struct):
            return
        else:
            self._dst_group_struct = np.asarray(dst_group_struct)
            self._dst_ngroups = len(dst_group_struct) - 1
            self._src_to_dst_matrix = partial_energy_matrix(
                dst_group_struct, self._src_group_struct
            )

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
            self.rxcache[rxkey] = (
                None if self.fullyloaded else self._load_reaction(nuc, rx, temp)
            )
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
        dst_sigma = (
            None
            if src_sigma is None
            else group_collapse(
                src_sigma, src_phi_g, dst_phi_g, self._src_to_dst_matrix
            )
        )
        return dst_sigma

    def shield_weights(self, num_dens, temp):
        """Builds the weights used during the self shielding calculations.
        Parameters
        ----------
        mat : array of floats.
            A map of the number densities of each of the nuclides of the material for
            which self-shielding is being calculated.
        data_source: pyne data_source
            Contains the cross section information for the isotopes in the material
            that is experiencing the self shielding.
        """
        reactions = {}
        for i in num_dens:
            rx = self.reaction(i, "total", temp)
            reactions[i] = 0.0 if rx is None else rx
        weights = {}
        for i in num_dens:
            weights[i] = 0.0
            for j in reactions:
                if j != i:
                    weights[i] += num_dens[j] * reactions[j]
            weights[i] = 1.0 / (weights[i] / num_dens[i] + reactions[i])
        self.slf_shld_wgts = weights

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
        """Loads the entire data source into memory."""
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
        return np.zeros(self.src_ngroups, dtype="f8")

    def discretize(self, nuc, rx, temp=300.0, src_phi_g=None, dst_phi_g=None):
        """Returns zeros."""
        return np.zeros(self.dst_ngroups, dtype="f8")

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

    _rx_avail = {
        rxname.id("total"): "t",
        rxname.id("scattering"): "s",
        rxname.id("elastic"): "e",
        rxname.id("inelastic"): "i",
        rxname.id("absorption"): "a",
        rxname.id("gamma"): "gamma",
        rxname.id("fission"): "f",
        rxname.id("alpha"): "alpha",
        rxname.id("proton"): "proton",
        rxname.id("deut"): "deut",
        rxname.id("trit"): "trit",
        rxname.id("z_2n"): "2n",
        rxname.id("z_3n"): "3n",
        rxname.id("z_4n"): "4n",
    }

    def __init__(self, **kwargs):
        super(SimpleDataSource, self).__init__(**kwargs)

    @property
    def exists(self):
        if self._exists is None:
            with tb.open_file(nuc_data, "r") as f:
                self._exists = "/neutron/simple_xs" in f
        return self._exists

    _USES_TEMP = False

    def _load_group_structure(self):
        """Sets the simple energy bounds array, E_g."""
        self.src_group_struct = np.array([14.0, 1.0, 2.53e-8, 0.0], dtype="float64")

    def _load_reaction(self, nuc, rx, temp=300.0):
        if rx not in self._rx_avail:
            return None
        cond = "nuc == {0}".format(nuc)
        sig = "sigma_" + self._rx_avail[rx]
        with tb.open_file(nuc_data, "r") as f:
            simple_xs = f.root.neutron.simple_xs
            fteen = [row[sig] for row in simple_xs.fourteen_MeV.where(cond)]
            fissn = [row[sig] for row in simple_xs.fission_spectrum_ave.where(cond)]
            therm = [row[sig] for row in simple_xs.thermal.where(cond)]
            if 0 == len(therm):
                therm = [row[sig] for row in simple_xs.thermal_maxwell_ave.where(cond)]
            if 0 == len(fteen) and 0 == len(fissn) and 0 == len(therm):
                rxdata = None
            else:
                rxdata = np.array([fteen[0], fissn[0], therm[0]], dtype="float64")
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
        dst_sigma = (src_sigma[2] * np.sqrt(2.53e-8)) / np.sqrt(center_g)
        dst_fissn = ((src_sigma[0] - src_sigma[1]) / 13.0) * (
            center_g - 1.0
        ) + src_sigma[1]
        mask = dst_sigma < dst_fissn
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
            self._dst_centers = (
                self._dst_group_struct[1:] + self._dst_group_struct[:-1]
            ) / 2.0
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
    _rx_avail = {
        rxname.id("np_1"): "np *",
        rxname.id("a_1"): "a  *",
        rxname.id("He3_1"): "h  *",
        rxname.id("z_2p_1"): "2p *",
        rxname.id("z_3n_1"): "3n *",
        rxname.id("d_1"): "d  *",
        rxname.id("npd"): "np/d",
        rxname.id("na"): "na",
        rxname.id("excited"): "*",
        rxname.id("nd"): "nd",
        rxname.id("gamma_1"): "g  *",
        rxname.id("z_3n"): "3n",
        rxname.id("np"): "np",
        rxname.id("nt"): "nt",
        rxname.id("t"): "t",
        rxname.id("nt_1"): "nt *",
        rxname.id("z_4n_1"): "4n *",
        rxname.id("na_1"): "na *",
        rxname.id("nd_1"): "nd *",
        rxname.id("t_1"): "t  *",
        rxname.id("a"): "a",
        rxname.id("z_2p"): "2p",
        rxname.id("d"): "d",
        rxname.id("gamma"): "g",
        rxname.id("He3"): "h",
        rxname.id("n"): "n",
        rxname.id("z_4n"): "4n",
        rxname.id("p"): "p",
        rxname.id("n_1"): "n  *",
        rxname.id("z_2a"): "2a",
        rxname.id("z_2n_1"): "2n *",
        rxname.id("z_2n"): "2n",
        rxname.id("nHe3_1"): "nh *",
        rxname.id("p_1"): "p  *",
        # not real or unique absorption reactions
        # rxname.id(''): "",
        # rxname.id(''): 'x',
        # rxname.id(''): 'x  *',
        # rxname.id(''): 'c',
        # rxname.id('fission'): 'f',
    }

    _USES_TEMP = False

    def __init__(self, **kwargs):
        super(CinderDataSource, self).__init__(**kwargs)

    def _load_group_structure(self):
        """Loads the cinder energy bounds array, E_g, from nuc_data."""
        with tb.open_file(nuc_data, "r") as f:
            E_g = np.array(f.root.neutron.cinder_xs.E_g)
        self.src_group_struct = E_g

    @property
    def exists(self):
        if self._exists is None:
            with tb.open_file(nuc_data, "r") as f:
                self._exists = "/neutron/cinder_xs" in f
        return self._exists

    def _load_reaction(self, nuc, rx, temp=300.0):
        fissrx = rxname.id("fission")
        absrx = rxname.id("absorption")

        # Set query condition
        if rx in self._rx_avail:
            if py3k:
                cond = "(from_nuc == {0}) & (reaction_type == {1})"
            else:
                cond = "(from_nuc == {0}) & (reaction_type == '{1}')"
            cond = cond.format(nuc, self._rx_avail[rx].encode())
        elif rx == fissrx:
            cond = "nuc == {0}".format(nuc)
        elif rx == absrx:
            cond = "(from_nuc == {0}) & (reaction_type != b'c')".format(nuc)
        else:
            return None

        # read & collapse data
        with tb.open_file(nuc_data, "r") as f:
            node = (
                f.root.neutron.cinder_xs.fission
                if rx == fissrx
                else f.root.neutron.cinder_xs.absorption
            )
            rows = [np.array(row["xs"]) for row in node.where(cond)]

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
    _rx_avail = {
        rxname.id("gamma"): b"1020",
        rxname.id("gamma_1"): b"1021",
        rxname.id("gamma_2"): b"1022",
        rxname.id("p"): b"1030",
        rxname.id("p_1"): b"1031",
        rxname.id("p_2"): b"1032",
        rxname.id("d"): b"1040",
        rxname.id("d_1"): b"1041",
        rxname.id("d_2"): b"1042",
        rxname.id("t"): b"1050",
        rxname.id("t_1"): b"1051",
        rxname.id("t_2"): b"1052",
        rxname.id("He3"): b"1060",
        rxname.id("He3_1"): b"1061",
        rxname.id("He3_2"): b"1062",
        rxname.id("a"): b"1070",
        rxname.id("a_1"): b"1071",
        rxname.id("a_2"): b"1072",
        rxname.id("z_2a"): b"1080",
        rxname.id("z_2p"): b"1110",
        rxname.id("z_2p_1"): b"1111",
        rxname.id("z_2p_2"): b"1112",
        rxname.id("z_2n"): b"160",
        rxname.id("z_2n_1"): b"161",
        rxname.id("z_2n_2"): b"162",
        rxname.id("z_3n"): b"170",
        rxname.id("z_3n_1"): b"171",
        rxname.id("z_3n_2"): b"172",
        rxname.id("fission"): b"180",
        rxname.id("na"): b"220",
        rxname.id("na_1"): b"221",
        rxname.id("na_2"): b"222",
        rxname.id("z_2na"): b"240",
        rxname.id("np"): b"280",
        rxname.id("np_1"): b"281",
        rxname.id("np_2"): b"282",
        rxname.id("n2a"): b"290",
        rxname.id("nd"): b"320",
        rxname.id("nd_1"): b"321",
        rxname.id("nd_2"): b"322",
        rxname.id("nt"): b"330",
        rxname.id("nt_1"): b"331",
        rxname.id("nt_2"): b"332",
        rxname.id("nHe3"): b"340",
        rxname.id("nHe3_1"): b"341",
        rxname.id("nHe3_2"): b"342",
        rxname.id("z_4n"): b"370",
        rxname.id("z_4n_1"): b"371",
        rxname.id("n"): b"40",
        rxname.id("n_1"): b"41",
        rxname.id("n_2"): b"42",
        rxname.id("z_3np"): b"420",
    }

    _avail_rx = dict([_[::-1] for _ in _rx_avail.items()])

    _USES_TEMP = False

    def __init__(self, **kwargs):
        super(EAFDataSource, self).__init__(**kwargs)

    def _load_group_structure(self):
        """Loads the EAF energy bounds array, E_g, from nuc_data."""
        with tb.open_file(nuc_data, "r") as f:
            E_g = np.array(f.root.neutron.eaf_xs.E_g)
        self.src_group_struct = E_g

    @property
    def exists(self):
        if self._exists is None:
            with tb.open_file(nuc_data, "r") as f:
                self._exists = "/neutron/eaf_xs" in f
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
        absrx = rxname.id("absorption")

        if rx in self._rx_avail:
            if py3k is True:
                cond = "(nuc_zz == {0}) & (rxnum == {1})"
            else:
                cond = "(nuc_zz == {0}) & (rxnum == '{1}')"
            cond = cond.format(nuc, self._rx_avail[rx])
        elif rx == absrx:
            cond = "(nuc_zz == {0})".format(nuc)
        else:
            return None

        # Grab data
        with tb.open_file(nuc_data, "r") as f:
            node = f.root.neutron.eaf_xs.eaf_xs
            rows = node.read_where(cond)
            # rows = [np.array(row['xs']) for row in node.where(cond)]

        if len(rows) == 0:
            rxdata = None
        elif 1 < len(rows):
            xss = rows["xs"]
            rxnums = rows["rxnum"]
            for rxnum, xs in zip(rxnums, xss):
                self.rxcache[nuc, self._avail_rx[rxnum]] = xs
            rxdata = xss.sum(axis=0)
        else:
            rxdata = rows[0]["xs"]

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
        absrx = rxname.id("absorption")
        with tb.open_file(nuc_data, "r") as f:
            node = f.root.neutron.eaf_xs.eaf_xs
            for row in node:
                nuc = row["nuc_zz"]
                rx = avail_rx[row["rxnum"]]
                xs = row["xs"]
                rxcache[nuc, rx] = xs
                abskey = (nuc, absrx)
                rxcache[abskey] = xs + rxcache.get(abskey, 0.0)
        self.fullyloaded = True


class ENDFDataSource(DataSource):
    """Evaluated Nuclear Data File cross section data source.  The ENDF file
    must exist for this data source to exist.

    Parameters
    ----------
    fh : string, file handle
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
            self.library = endf.Library(fh)
        self.rxcache = {}
        self.dst_group_struct = dst_group_struct
        self._src_phi_g = src_phi_g

    @property
    def exists(self):
        if self._exists is None:
            if isinstance(self.fh, basestring):
                self._exists = os.path.isfile(self.fh)
            else:
                self._exists = isinstance(self.fh, IO_TYPES)
        return self._exists

    def _load_group_structure(self, nuc, rx, nuc_i=None):
        """Loads the group structure from ENDF file."""
        self.library._read_res(nuc)
        mt = rxname.mt(rx)
        rxdata = self.rxcache[nuc, rx, nuc_i]
        xsdata = self.library.get_xs(nuc, rx, nuc_i)[0]
        intpoints = xsdata["intpoints"]  # [::-1]
        e_int = xsdata["e_int"]
        src_group_struct = [e_int[intpoint - 1] for intpoint in intpoints[::-1]]
        rxdata["src_group_struct"] = src_group_struct
        rxdata["src_ngroups"] = len(intpoints)
        rxdata["src_phi_g"] = (
            np.ones(len(intpoints), dtype="f8")
            if rxdata["_src_phi_g"] is None
            else np.asarray(rxdata["src_phi_g"])
        )

    def _load_reaction(self, nuc, rx, nuc_i, src_phi_g=None, temp=300.0):
        """Note: ENDF data does not use temperature information (temp)

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
        rx = rxname.mt(rx)
        if nuc_i not in self.library.structure[nuc]["data"]:
            self.library._read_res(nuc)
        rxdict = self.library.get_xs(nuc, rx, nuc_i)[0]
        rxdict["_src_phi_g"] = src_phi_g
        self.rxcache[nuc, rx, nuc_i] = rxdict
        self._load_group_structure(nuc, rx, nuc_i)
        return rxdict

    def reaction(self, nuc, rx, nuc_i=None):
        """Get reaction data.

        Parameters
        ----------
        nuc : int or str
            Nuclide containing the reaction.
        rx : int or str
            Desired reaction
        nuc_i : int or str
            Nuclide containing the reaction. Defaults to nuc.
        group_bounds : tuple
            Low and high energy bounds of the new group.

        Returns
        ----------
        rxdict : dict
            Dictionary containing source group structure, energy values, cross-
            section data, and interpolation data."""
        if nuc_i is None:
            nuc_i = nuc
        nuc = nucname.id(nuc)
        rx = rxname.mt(rx)
        nuc_i = nucname.id(nuc_i)
        if (nuc, rx, nuc_i) not in self.rxcache:
            self._load_reaction(nuc, rx, nuc_i)
        return self.rxcache[nuc, rx, nuc_i]

    def discretize(self, nuc, rx, temp=300.0, src_phi_g=None, dst_phi_g=None):
        """Discretizes the reaction channel from the source group structure to that
        of the destination weighted by the group fluxes.

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
        nuc = nucname.id(nuc)
        nuc_i = nucname.id(nuc)
        rx = rxname.mt(rx)
        rxdata = self.reaction(nuc, rx, nuc_i=nuc_i)
        xs = rxdata["xs"]
        dst_group_struct = rxdata["dst_group_struct"]
        intpoints = [intpt for intpt in rxdata["intpoints"][::-1]]
        intschemes = rxdata["intschemes"][::-1]
        e_int = rxdata["e_int"]

        src_bounds = [e_int[intpoint - 1] for intpoint in intpoints]
        src_dict = dict(zip(src_bounds, intschemes))
        dst_bounds = zip(dst_group_struct[1:], dst_group_struct[:-1])
        dst_sigma = [
            self.integrate_dst_group(dst_bound, src_bounds, src_dict, e_int, xs)
            for dst_bound in dst_bounds
        ]
        return dst_sigma

    def integrate_dst_group(self, dst_bounds, src_bounds, src_dict, e_int, xs):
        dst_low, dst_high = dst_bounds
        src_bounds = np.array(src_bounds)

        # We're going to have to integrate over each zone bounded by the edges
        # of a destination bin or by the edges of a source bin
        internal_src_bounds = [bd for bd in src_bounds if dst_low < bd < dst_high]
        integration_bounds = [dst_low, dst_high]
        integration_bounds[1:1] = internal_src_bounds
        integration_bounds = list(zip(integration_bounds[:-1], integration_bounds[1:]))

        schemes = [
            src_dict[src_bounds[src_bounds >= high][0]]
            for low, high in integration_bounds
        ]

        integration_args = [
            (scheme, e_int, xs, bds[0], bds[1])
            for scheme, bds in zip(schemes, integration_bounds)
        ]
        return sum(
            [self._integrate_range_nonnormal(*args) for args in integration_args]
        ) / (dst_high - dst_low)

    def _integrate_range_nonnormal(self, scheme, e_int, xs, low, high):
        """De-normalizes the integral over a certain range. Useful when the
        range is integrated piecewise over several chunks - you don't want
        each chunk to be individually normalized.

        Parameters
        ----------
        scheme : int
            ENDF-coded interpolation scheme to be used between the data points.
        e_int : ndarray
            Array of energy values to integrate over.
        xs : ndarray
            Array of cross-sections corresponding to e_int.
        low, high : float
            Lower and upper bounds of integration.

        Returns
        -------
        sigma * dE : float
            Non-normalized integral.

        """
        dE = high - low
        sigma = self.library.integrate_tab_range(scheme, e_int, xs, low, high)
        return sigma * dE


class OpenMCDataSource(DataSource):
    """Data source for ACE data that is listed in an OpenMC cross_sections.xml
    file. This data source discretizes the reactions to a given group
    stucture when the reactions are loaded in. Reseting this source group
    structure will clear the reaction cache.
    """

    self_shield_reactions = {
        rxname.id("fission"),
        rxname.id("gamma"),
        rxname.id("total"),
    }

    def __init__(self, cross_sections=None, src_group_struct=None, **kwargs):
        """Parameters
        ----------
        cross_sections : openmc_utils.CrossSections or string or file-like, optional
            Path or file to OpenMC cross_sections.xml
        src_group_struct : array-like, optional
            The group structure to discretize the ACE data to, defaults to
            ``np.logspace(1, -9, 101)``.
        kwargs : optional
            Keyword arguments to be sent to DataSource base class.

        """
        if not isinstance(cross_sections, openmc_utils.CrossSections):
            cross_sections = cross_sections or os.getenv("CROSS_SECTIONS")
            cross_sections = openmc_utils.CrossSections(f=cross_sections)
        self.cross_sections = cross_sections
        self._src_group_struct = src_group_struct
        super(OpenMCDataSource, self).__init__(**kwargs)
        self.libs = {}  # cross section libraries, index by openmc_utils.AceTables

    @property
    def exists(self):
        if self._exists is None:
            self._exists = len(self.cross_sections.ace_tables) > 0
        return self._exists

    def _load_group_structure(self):
        if self._src_group_struct is None:
            self._src_group_struct = np.logspace(1, -9, 101)
        self.src_group_struct = self._src_group_struct

    def pointwise(self, nuc, rx, temp=300.0):
        """Returns pointwise reaction data from ACE files indexed by OpenMC.

        Parameters
        ----------
        nuc : int
            Nuclide id.
        rx : int
            Reaction id.
        temp : float, optional
            The nuclide temperature in [K].

        Returns
        -------
        E_points : array-like
            The array or energy points that the reaction is evaluated at.
        rawdata : array-like
            Raw pointwise reaction data.
        """
        nuc = nucname.id(nuc)
        rx = rxname.id(rx)
        try:
            mt = rxname.mt(rx)
        except RuntimeError:
            return None
        totrx = rxname.id("total")
        absrx = rxname.id("absorption")
        ace_tables = self._rank_ace_tables(nuc, temp=temp)
        lib = ntab = None
        for atab in ace_tables:
            if os.path.isfile(atab.abspath or atab.path):
                if atab not in self.libs:
                    lib = self.libs[atab] = ace.Library(atab.abspath or atab.path)
                    lib.read(atab.name)
                lib = self.libs[atab]
                ntab = lib.tables[atab.name]
                if mt in ntab.reactions or rx == totrx or rx == absrx:
                    break
                lib = ntab = None
        if lib is None:
            return None  # no reaction available
        E_g = self.src_group_struct
        E_points = ntab.energy
        if rx == totrx:
            rawdata = ntab.sigma_t
        elif rx == absrx:
            rawdata = ntab.sigma_a
        else:
            ntabrx = ntab.reactions[mt]
            if ntabrx.IE is None or ntabrx.IE == 0:
                rawdata = ntabrx.sigma
            else:
                rawdata = np.empty(len(E_points), dtype="f8")
                rawdata[: ntabrx.IE] = 0.0
                rawdata[ntabrx.IE :] = ntabrx.sigma
        if (E_g[0] <= E_g[-1] and E_points[-1] <= E_points[0]) or (
            E_g[0] >= E_g[-1] and E_points[-1] >= E_points[0]
        ):
            E_points = E_points[::-1]
            rawdata = rawdata[::-1]
        return E_points, rawdata

    def _load_reaction(self, nuc, rx, temp=300.0):
        """Loads reaction data from ACE files indexed by OpenMC.

        Parameters
        ----------
        nuc : int
            Nuclide id.
        rx : int
            Reaction id.
        temp : float, optional
            The nuclide temperature in [K].
        """
        rtn = self.pointwise(nuc, rx, temp=temp)
        if rtn is None:
            return
        E_points, rawdata = rtn
        E_g = self.src_group_struct
        if self.atom_dens.get(nuc, 0.0) > 1.0e19 and rx in self.self_shield_reactions:
            rxdata = self.self_shield(nuc, rx, temp, E_points, rawdata)
        else:
            rxdata = bins.pointwise_linear_collapse(E_g, E_points, rawdata)
        return rxdata

    def self_shield(self, nuc, rx, temp, E_points, xs_points):
        """Calculates the self shielded cross section for a given nuclide
        and reaction. This calculation uses the Bonderanko method.

        Parameters
        ----------
        nuc : int
            Nuclide id.
        rx : int
            Reaction id.
        temp : float, optional
            The nuclide temperature in [K].
        E_points : array like
            The point wise energies.
        xs_points : array like
            Point wise cross sections

        Returns
        -------
        rxdata : array like
            collapsed self shielded cross section for nuclide nuc and reaction
            rx
        """
        sigb = self.bkg_xs(nuc, temp=temp)
        e_n = self.src_group_struct
        sig_b = np.ones(len(E_points), "f8")
        for n in range(len(sigb)):
            sig_b[(e_n[n] <= E_points) & (E_points <= e_n[n + 1])] = sigb[n]
        rtn = self.pointwise(nuc, "total", temp)
        if rtn is None:
            sig_t = 0.0
        else:
            sig_t = rtn[1]
        numer = bins.pointwise_linear_collapse(
            self.src_group_struct, E_points, xs_points / (E_points * (sig_b + sig_t))
        )
        denom = bins.pointwise_linear_collapse(
            self.src_group_struct, E_points, 1.0 / (E_points * (sig_b + sig_t))
        )
        return numer / denom

    def bkg_xs(self, nuc, temp=300):
        """Calculates the background cross section for a nuclide (nuc)

        Parameters
        ----------
        nuc : int
            Nuclide id.
        temp : float, optional
            The nuclide temperature in [K].

        Returns
        -------
        sig_b : array like
            Group wise background cross sections.
        """
        e_n = self.src_group_struct
        sig_b = np.zeros(self.src_ngroups, float)
        for i, a in self.atom_dens.items():
            if i == nuc:
                continue
            rtn = self.pointwise(i, "total", temp)
            if rtn is None:
                continue
            sig_b += a * bins.pointwise_linear_collapse(e_n, rtn[0], rtn[1])
        return sig_b / self.atom_dens.get(nuc, 0.0)

    def _rank_ace_tables(self, nuc, temp=300.0):
        """Filters and sorts the potential ACE tables based on nucliude and
        temperature.
        """
        tabs = [t for t in self.cross_sections.ace_tables if t.nucid == nuc]
        if len(tabs) == 0:
            return tabs
        temps = {t.temperature for t in tabs}
        temps = sorted(temps, key=lambda s: abs(float(s) - temp * MeV_per_K))
        nearest_temp = temps[0]
        tabs = [t for t in tabs if t.temperature == nearest_temp]
        tabs.sort(reverse=True, key=lambda t: t.name)
        return tabs

    def load(self, temp=300.0):
        """Loads the entire data source into memory. This can be expensive for
        lots of ACE data.

        Parameters
        ----------
        temp : float, optional
            Temperature [K] of material, defaults to 300.0.

        """
        for atab in self.cross_sections.ace_tables:
            if os.path.isfile(atab.abspath or atab.path):
                lib = self.libs[atab] = ace.Library(atab.abspath or atab.path)
                lib.read(atab.name)


class StatePointDataSource(DataSource):
    """Data source for reactions coming from OpenMC state points"""

    def __init__(self, state_point, tallies, num_den, phi_tot, **kwargs):
        """Parameters
        ----------
        state_point : string
            Path to the OpenMC statepoint file to be used to build the data_source
        tallies : array-like
            The tally id's used to pull the cross sections from.
        num_dens: map of int to float
            A map containing the number densities of the nuclides in the
            material used in the statepoint.
        phi_tot: array of floats
            The total flux within the reactor
        kwargs : optional
            Keyword arguments to be sent to DataSource base class.

        """
        self.state_point = state_point
        self.tallies = tallies
        self.particles = state_point.n_particles
        self.reactions = {}
        self._load_reactions(num_den, phi_tot)
        super(StatePointDataSource, self).__init__(**kwargs)

    @property
    def exists(self):
        return True

    def _load_group_structure(self):
        """Loads the group structure from a tally in OpenMC. It is
        assumed that all tallies have the same group structure
        """
        self._src_group_struct = (
            self.state_point.tallies[self.tallies[0]].filters[0].bins[::-1]
        )
        self.src_group_struct = self._src_group_struct

    def _load_reactions(self, num_dens, phi_tot):
        """Loads the group structure from a tally in OpenMC. It is
        assumed that all tallies have the same group structure

        Parameters
        ----------
        state_point: OpenMC statepoint file
            The statepoint file that will be used to load the reaction
            rates to determine the microscopic cross sections for the
            statepoint.
        num_dens: map of int to float
            A map containing the number densities of the nuclides in the
            material used in the statepoint.
        phi_tot: array of floats
            The total flux within the reactor
        """
        for tally in self.tallies:
            sp = self.state_point.tallies[tally]
            for nuclide in sp.nuclides:
                for score in sp.scores:
                    self.reactions[nucname.id(nuclide.name), score] = (
                        sp.get_values([score], [], [], [nuclide.name]).flatten()
                        * self.particles
                    )
                    weight = 1.0e24 / abs(phi_tot * num_dens[nucname.id(nuclide.name)])
                    self.reactions[nucname.id(nuclide.name), score] *= weight

    def reaction(self, nuc, rx, temp):
        """Loads reaction data from ACE files indexed by OpenMC.

        Parameters
        ----------
        nuc : int
            Nuclide id.
        rx : int
            Reaction id.
        Return
        ------
        Array containing the cross sections for the reaction
        requested.
        """
        rxkey = (nuc, rx)
        if rxkey not in self.reactions:
            return None
        else:
            return self.reactions[rxkey]
