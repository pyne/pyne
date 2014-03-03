"""This module provides a generic Transmuter interface.
"""

import numpy as np

class Transmuter(object):
    def __init__(self, t=0.0, phi=0.0, temp=300.0, tol=1e-7, xscache=None):
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
        xscache : XSCache, optional
            A cross section cache to generate cross sections with.
        """
        self.xscache = xscache
        self.t = t
        self._phi = None
        self.phi = phi
        self.temp = temp
        self.tol = tol

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
        if self.xscache is not None:
            for ds in self.xscache.data_sources:
                ds.src_phi_g = flux
            self.xscache['phi_g'] = np.array([flux.sum()])
        self._phi = flux

    def transmute(*args, **kwargs):
        """All derived classes must implement a transmute method."""
        pass

    def transmute_mesh(self, *args, **kwargs):
        """All derived classes must implement a transmute_mesh method."""
        pass
