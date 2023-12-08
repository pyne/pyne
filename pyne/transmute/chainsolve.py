"""This module implements an ALARA[1]-like chain-based transmutation solver.

   [1] Wilson, P. P. H. "ALARA: Analytic Laplacian Adaptive Radioactivity 
   Analysis," a Ph.D. Dissertation, University of Wisconsin, Madison, WI, 1999.
"""
from __future__ import division
from pyne.utils import QA_warn

import numpy as np
from scipy import linalg

# from scipy import sparse  # <-- SPARSE

from pyne import utils
from pyne import data
from pyne import rxname
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material, from_atom_frac
from pyne.xs.data_source import NullDataSource, EAFDataSource
from pyne.xs.cache import XSCache
from pyne.xs.channels import sigma_a

QA_warn(__name__)


class Transmuter(object):
    """A class for transmuting materials using an ALARA-like chain solver."""

    def __init__(
        self, t=0.0, phi=0.0, temp=300.0, tol=1e-7, rxs=None, log=None, *args, **kwargs
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
        rxs : set of ints or strs
            Reaction ids or names to use in transmutation that produce well-defined
            children.  This set should thus not include fission.  If None, then the
            reactions from EAF are used.
        log : file-like or None
            The log file object should be written. A None imples the log is
            not desired.
        args : tuple, optional
            Other arguments ignored for compatibility with other Transmuters.
        kwargs : dict, optional
            Other keyword arguments ignored for compatibility with other Transmuters.
        """
        eafds = EAFDataSource()
        eafds.load(temp=temp)
        gs = np.array([eafds.src_group_struct[0], eafds.src_group_struct[-1]])
        eafds.dst_group_struct = gs
        self.xscache = XSCache(
            group_struct=gs,
            data_sources=(
                eafds,
                NullDataSource,
            ),
        )

        self.t = t
        self._phi = None
        self.phi = phi
        self.temp = temp
        self.log = log
        self.tol = tol
        if rxs is None:
            rxs = [
                "gamma",
                "gamma_1",
                "gamma_2",
                "p",
                "p_1",
                "p_2",
                "d",
                "d_1",
                "d_2",
                "t",
                "t_1",
                "t_2",
                "He3",
                "He3_1",
                "He3_2",
                "a",
                "a_1",
                "a_2",
                "z_2a",
                "z_2p",
                "z_2p_1",
                "z_2p_2",
                "z_2n",
                "z_2n_1",
                "z_2n_2",
                "z_3n",
                "z_3n_1",
                "z_3n_2",
                "na",
                "na_1",
                "na_2",
                "z_2na",
                "np",
                "np_1",
                "np_2",
                "n2a",
                "nd",
                "nd_1",
                "nd_2",
                "nt",
                "nt_1",
                "nt_2",
                "nHe3",
                "nHe3_1",
                "nHe3_2",
                "z_4n",
                "z_4n_1",
                "n",
                "n_1",
                "n_2",
                "z_3np",
            ]
        rxs = set([rxname.id(rx) for rx in rxs])
        rxs.discard(rxname.id("fission"))
        self.rxs = rxs

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

    def transmute(self, x, t=None, phi=None, tol=None, log=None, *args, **kwargs):
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
        mw_x = x.molecular_mass()
        y = from_atom_frac(y_atoms, atoms_per_molecule=x.atoms_per_molecule)
        # even though it doesn't look like it, the following line is actually
        #   mass_y = MW_y * mass_x / MW_x
        y.mass *= x.mass / mw_x
        return y

    def _transmute_partial(self, nuc):
        """Core method to transmute a material into its daughters.
        This method assumes that the initial nuclide has unit density.

        Parameters
        ----------
        nuc : int
            Nuclide id to be transmuted.

        Returns
        -------
        partial : dict
            A dictionary containing number densities for each nuclide after
            the transmutation is carried out for the input nuclide. Keys are
            nuclide ids and values are float number densities for the coupled
            # what is the coupled?.
        """
        dest = self._get_destruction(nuc)
        # DENSE
        A = np.empty((1, 1), float)
        A[0, 0] = -dest
        # A = sparse.csr_matrix([[-dest]])  # <-- SPARSE
        rootval = np.exp(-dest * self.t)
        partial = {nuc: rootval}
        self._traversal(nuc, A, partial)
        return partial

    def _get_destruction(self, nuc, decay=True):
        """Computes the destruction rate of the nuclide.

        Parameters
        ----------
        nuc : int
            Name of the nuclide in question
        decay : bool
            True if the decay constant should be added to the returned value.
            False if only destruction from neutron reactions should be considered.

        Returns
        -------
        d : float
            Destruction rate of the nuclide.

        """
        xscache = self.xscache
        sig_a = sigma_a(nuc, xs_cache=xscache)
        d = utils.from_barns(sig_a[0], "cm2") * xscache["phi_g"][0]
        if decay and not np.isnan(data.decay_const(nuc)):
            d += data.decay_const(nuc)
        return d

    def _grow_matrix(self, A, prod, dest):
        """Grows the given matrix by one row and one column, adding necessary
        production and destruction rates.

        Parameters
        ----------
        A : NumPy 2-dimensional array
            The original matrix that must be grown.
        prod : float
            The production rate of the next nuclide in the chain.
        dest : float
            The destruction rate of the next nuclide in the chain.

        Returns
        -------
        B : NumPy 2-dimensional array
            The grown matrix
        """
        shape = A.shape
        n = shape[0]
        # DENSE
        # Add row and column to current matrix
        B = np.empty((n + 1, n + 1), dtype=float)
        B[:n, :n] = A
        B[n, : n - 1] = 0
        B[:n, n] = 0
        # Update new matrix with provided data
        B[n, n - 1] = prod
        B[n, n] = -dest
        # SPARSE
        # B = sparse.bmat([[A, None], [None, [[-dest]]]]).tocsr()
        # B[n,n-1] = prod
        return B

    def _traversal(self, nuc, A, out, depth=0):
        """Nuclide transmutation traversal method.

        This method will traverse the reaction tree recursively, using a DFS
        algorithm. On termination, the method will return all number densities
        after a given time that are a result of the starting nuclide.

        Parameters
        ----------
        nuc : int
            ID of the active nuclide for the traversal.
        A : NumPy 2-dimensional array
            Current state of the coupled equation matrix.
        out : dict
            A dictionary containing the final recorded number densities for each
            nuclide. Keys are nuclide names in integer id form. Values are
            number densities for the coupled nuclide in float format.  This is
            modified in place.
        depth : int
            Current depth of traversal (root at 0). Should never be provided by user.

        """
        t = self.t
        tol = self.tol
        phi = self.xscache["phi_g"][0]
        temp = self.temp
        xscache = self.xscache
        if self.log is not None:
            self._log_tree(depth, nuc, 1.0)
        prod = {}
        # decay info
        lam = data.decay_const(nuc)
        decay_branches = {} if lam == 0 else self._decay_branches(nuc)
        for decay_child, branch_ratio in decay_branches.items():
            prod[decay_child] = lam * branch_ratio
        # reaction daughters
        for rx in self.rxs:
            try:
                child = rxname.child(nuc, rx)
            except RuntimeError:
                continue
            child_xs = xscache[nuc, rx, temp][0]
            rr = utils.from_barns(child_xs, "cm2") * phi  # reaction rate
            prod[child] = rr + prod.get(child, 0.0)
        # Cycle production dictionary
        for child in prod:
            # Grow matrix
            d = self._get_destruction(child)
            B = self._grow_matrix(A, prod[child], d)
            # Create initial density vector
            n = B.shape[0]
            N0 = np.zeros((n, 1), dtype=float)
            N0[0] = 1.0
            # Compute matrix exponential and dot with density vector
            eB = linalg.expm(B * t)
            N_final = np.dot(eB, N0)  # <-- DENSE
            # N_final = eB.dot(N0)  # <-- SPARSE
            if self.log is not None:
                self._log_tree(depth + 1, child, N_final[-1])
            # Check against tolerance and continue traversal
            if N_final[-1] > tol:
                self._traversal(child, B, out, depth=depth + 1)
            # On recursion exit or truncation, write data from this nuclide
            outval = N_final[-1, 0] + out.get(child, 0.0)
            if 0.0 < outval:
                out[child] = outval

    def _log_tree(self, depth, nuc, numdens):
        """Logging method to track path of _traversal.

        Parameters
        ----------
        depth : integer
            Current depth of traversal (root at 0).
        nuc : nucname
            Current nuclide in traversal.
        numdens : float
            Current number density of nuc.
        tree : File
            File to write tree log to.

        """
        # Don't print a zero density.
        if numdens == 0.0:
            return
        space = "   |"
        entry = "{spacing}--> {name} {numdens}\n".format(
            spacing=depth * space, numdens=numdens, name=nucname.name(nuc)
        )
        self.log.write(entry)
        self.log.flush()

    def _decay_branches(self, nuc):
        """Returns a dictionary that contains the decay children of nuc as keys
        to the branch ratio of that child's decay process.

        Parameters
        ----------
        nuc : int
            Name of parent nuclide to get decay children of.

        Returns
        -------
        decay_branches : dictionary
            Keys are decay children of nuc in zzaaam format.
            Values are the branch ratio of the decay child.

        """
        decay_branches = {}
        children = data.decay_children(nuc)
        for child in children:
            decay_branches[child] = data.branch_ratio(nuc, child)
        return decay_branches
