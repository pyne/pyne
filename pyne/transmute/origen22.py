"""This module implements an ORIGEN v2.2 transmutation solver.
"""
from __future__ import print_function, division

import os
import subprocess
import tempfile

try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping
from pyne.utils import QA_warn

import numpy as np

from pyne import data
from pyne import rxname
from pyne import nucname
from pyne import nuc_data
from pyne import origen22
from pyne import utils
from pyne.material import Material, from_atom_frac
from pyne.xs.data_source import NullDataSource, SimpleDataSource, EAFDataSource
from pyne.xs.cache import XSCache

QA_warn(__name__)


class Transmuter(object):
    """A class for transmuting materials using ORIGEN v2.2."""

    def __init__(
        self,
        t=0.0,
        phi=0.0,
        temp=300.0,
        tol=1e-7,
        cwd="",
        base_tape9=origen22.BASE_TAPE9,
        xscache=None,
        o2exe="o2_therm_linux.exe",
        *args,
        **kwargs
    ):
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
        cwd : str, optional
            Current working directory for origen runs. Defaults to this dir.
        base_tape9 : str or dict, optional
            A base TAPE9.INP file.  If this is a str it is interpreted as a path
            to a file, which is then read in and parsed.  If this is a dict, it is
            assumed to be in the format described in the main origen22 module.
        xscache : XSCache, optional
            A cross section cache to generate cross sections with.
        o2exe : str, optional
            Name or path to ORIGEN 2.2 executable.
        args : tuple, optional
            Other arguments ignored for compatibility with other Transmuters.
        kwargs : dict, optional
            Other keyword arguments ignored for compatibility with other Transmuters.
        """
        if not isinstance(base_tape9, Mapping):
            base_tape9 = origen22.parse_tape9(tape9=base_tape9)
        self.base_tape9 = base_tape9

        if xscache is None:
            eafds = EAFDataSource()
            eafds.load(temp=temp)
            gs = np.array([eafds.src_group_struct[0], eafds.src_group_struct[-1]])
            eafds.dst_group_struct = gs
            xscache = XSCache(
                group_struct=gs, data_source_classes=[SimpleDataSource, NullDataSource]
            )
            xscache.load(temp=temp)
            xscache.data_sources.insert(0, eafds)
        self.xscache = xscache

        self.t = t
        self._phi = None
        self.phi = phi
        self.temp = temp
        self.tol = tol
        self.cwd = os.path.abspath(cwd)
        self.o2exe = o2exe

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
        for ds in self.xscache.data_sources:
            ds.src_phi_g = flux
        self.xscache["phi_g"] = np.array([flux.sum()])
        self._phi = flux

    def transmute(
        self,
        x,
        t=None,
        phi=None,
        tol=None,
        cwd=None,
        xscache=None,
        o2exe=None,
        *args,
        **kwargs
    ):
        """Transmutes a material into its daughters.

        Parameters
        ----------
        x : Material or similar
            Input material for transmutation.
        t : float
            Transmutations time [sec].
        phi : float or array of floats
            Neutron flux vector [n/cm^2/sec].  Currently this must either be
            a scalar or match the group structure of EAF.
        tol : float
            Tolerance level for chain truncation.
        cwd : str, optional
            Current working directory for origen runs. Defaults to this dir.
        xscache : XSCache, optional
            A cross section cache to generate cross sections with.
        o2exe : str, optional
            Name or path to ORIGEN 2.2 executable.

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
        if tol is not None:
            self.tol = tol
        if cwd is not None:
            self.cwd = os.path.abspath(cwd)
        if xscache is not None:
            self.xscache = xscache
        if o2exe is not None:
            self.o2exe = o2exe

        # prepare new tape9
        nucs = set(x.comp.keys())
        base_tape9 = self.base_tape9
        decay_nlb, xsfpy_nlb = origen22.nlbs(base_tape9)
        new_tape9 = origen22.xslibs(nucs=nucs, xscache=self.xscache, nlb=xsfpy_nlb)
        t9 = origen22.merge_tape9([new_tape9, base_tape9])

        # write out files
        origen22.write_tape4(x, outfile=os.path.join(self.cwd, "TAPE4.INP"))
        origen22.write_tape5_irradiation(
            "IRF",
            self.t / 86400.0,
            self.xscache["phi_g"][0],
            outfile=os.path.join(self.cwd, "TAPE5.INP"),
            decay_nlb=decay_nlb,
            xsfpy_nlb=xsfpy_nlb,
            cut_off=self.tol,
        )
        origen22.write_tape9(t9, outfile=os.path.join(self.cwd, "TAPE9.INP"))

        # run origen & get results
        f = tempfile.NamedTemporaryFile()
        try:
            subprocess.check_call([self.o2exe], cwd=self.cwd, stdout=f, stderr=f)
        except subprocess.CalledProcessError:
            f.seek(0)
            print("ORIGEN output:\n\n{0}".format(f.read()))
            raise
        finally:
            f.close()
        t6 = origen22.parse_tape6(tape6=os.path.join(self.cwd, "TAPE6.OUT"))
        y = t6["materials"][-1]
        return y
