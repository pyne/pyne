"""This module implements an ALARA-like chain-based transmutation solver.
"""
from __future__ import print_function
import subprocess
from collections import Mapping

import numpy as np

from pyne import utils
from pyne import data
from pyne import rxname
from pyne import nucname
from pyne import nuc_data
from pyne import origen22
from pyne.material import Material, from_atom_frac
from pyne.xs.data_source import NullDataSource, EAFDataSource
from pyne.xs.cache import XSCache

class Transmuter(object):
    """A class for transmuting materials using an ALARA-like chain solver."""

    def __init__(self, t=0.0, phi=0.0, temp=300.0, tol=1e-7, 
                 base_tape9=origen22.BASE_TAPE9, *args, **kwargs):
        """Parameters
        ----------
        t : float
            Transmutations time [sec].
        phi : float or array of floats
            Neutron flux vector [n/cm^2/sec].  Currently this must either be 
            a scalar or match the group structure of EAF.
        temp : float, optional
            Temperature [K] of material, defaults to 300.0.
        tol : float
            Tolerance level for chain truncation.
        base_tape9 : str or dict, optional
            A base TAPE9.INP file.  If this is a str it is interpreted as a path 
            to a file, which is then read in and parsed.  If this is a dict, it is
            assumed to be in the format described in the main origen22 module.
        args : tuple, optional
            Other arguments ignored for compatibility with other Transmuters.
        kwargs : dict, optional
            Other keyword arguments ignored for compatibility with other Transmuters.
        """
        eafds = EAFDataSource()
        eafds.load(temp=temp)
        gs = np.array([eafds.src_group_struct[0], eafds.src_group_struct[-1]])
        eafds.dst_group_struct = gs
        self.xs_cache = XSCache(group_struct=gs, data_source_classes=(NullDataSource,))
        self.xs_cache.data_sources.insert(0, eafds)

        self.t = t
        self._phi = None
        self.phi = phi
        self.temp = temp
        self.tol = tol

        if not isinstance(base_tape9, Mapping):
            base_tape9 = origen22.parse_tape9(tape9=base_tape9)
        self.base_tape9 = base_tape9

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, flux):
        """Ensures that the flux is correctly formatted."""
        flux = np.asarray(flux)
        if flux.ndim == 0:
            _ = np.empty(175, float)
            _.fill(flux / 175.0)
            flux = _
        elif flux.ndim == 1 and flux.shape[0] != 175:
            raise ValueError("Group structure must match EAF.")
        elif flux.ndim > 1:
            raise ValueError("The flux vector must be 0- or 1-dimensional.")
        if not np.all(flux >= 0.0):
            raise ValueError("Flux entries must be non-negative.")
        for ds in self.xs_cache.data_sources:
            ds.src_phi_g = flux
        self.xs_cache['phi_g'] = np.array([flux.sum()])
        self._phi = flux

    def transmute(self, x, t=None, phi=None, tol=None, log=None):
        """Transmutes a material into its daughters.

        Parameters
        ----------
        x : Marterial or similar
            Input material for transmutation.
        t : float
            Transmutations time [sec].
        phi : float or array of floats
            Neutron flux vector [n/cm^2/sec].  Currently this must either be 
            a scalar or match the group structure of EAF.
        tol : float
            Tolerance level for chain truncation.
        log : file-like or None
            The log file object should be written. A None imples the log is 
            not desired.

        Returns
        -------
        y : Material
            The output material post-transmutation.

        """
        if not isinstance(x, Material):
            x = Material(x)
        if t is not None:
            self.t = t
        if phi is not None:
            self.phi = phi
        if log is not None:
            self.log = log
        if tol is not None:
            self.tol = tol

        x_atoms = x.to_atom_frac()
        y_atoms = {}
        for nuc, adens in x_atoms.items():
            # Find output for root of unit density and scale all output by 
            # actual nuclide density and add to final output.
            partial = self._transmute_partial(nuc)
            for part_nuc, part_adens in partial.items():
                y_atoms[part_nuc] = part_adens * adens + y_atoms.get(part_nuc, 0.0)
        mw_x = x.molecular_weight()
        y = from_atom_frac(y_atoms, atoms_per_mol=x.atoms_per_mol)
        # even though it doesn't look likt it, the following line is actually
        #   mass_y = MW_y * mass_x / MW_x
        y.mass *= x.mass / mw_x 
        return y

