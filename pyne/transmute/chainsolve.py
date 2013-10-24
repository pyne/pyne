"""This module implements an ALARA-like chain-based transmutation solver.
"""

import numpy as np
from scipy import linalg

from pyne import utils
from pyne import data
from pyne import rxname
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material, from_atom_frac
from pyne.xs.data_source import NullDataSource, EAFDataSource
from pyne.xs.cache import XSCache
from pyne.xs.channels import sigma_a

class Transmuter(object):
    """A class for transmuting materials using an ALARA-like chain solver."""

    def __init__(self, t=0.0, phi=0.0, temp=300.0, tol=1e-7, rxs=None, log=None):
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
            Reaction ids or names to use in transmutation which produce well-defined 
            children.  This set should thus not include fission.  If None, then the 
            reactions from EAF are used.
        log : file-like or None
            The log file object should be written. A None imples the log is 
            not desired.
        """
        eafds = EAFDataSource()
        gs = np.array([eafds.src_group_struct[0], eafds.src_group_struct[-1]])
        eafds.dst_group_struct = gs
        self.xs_cache = XSCache(group_struct=gs, data_source_classes=(NullDataSource))
        self.xs_cache.data_sources.insert(0, eafds)

        self.t = t
        self._phi = None
        self.phi = phi
        self.temp = temp
        self.log = log
        self.tol = tol
        if rxs is None:
            rxs = ['gamma', 'gamma_1', 'gamma_2', 'p', 'p_1', 'p_2', 'd', 'd_1', 
                   'd_2', 't', 't_1', 't_2', 'He3', 'He3_1', 'He3_2', 'a', 'a_1', 
                   'a_2', 'z_2a', 'z_2p', 'z_2p_1', 'z_2p_2', 'z_2n', 'z_2n_1', 
                   'z_2n_2', 'z_3n', 'z_3n_1', 'z_3n_2', 'na', 'na_1', 'na_2', 
                   'z_2na', 'np', 'np_1', 'np_2', 'n2a', 'nd', 'nd_1', 'nd_2', 
                   'nt', 'nt_1', 'nt_2', 'nHe3', 'nHe3_1', 'nHe3_2','z_4n', 
                   'z_4n_1', 'n', 'n_1', 'n_2', 'z_3np']
        rxs = set([rxname.id(rx) for rx in rxs])
        rxs.discard(rxname.id('fission'))
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
            if flux == 0.0
                _.fill(0.0)
            else:
                _.fill(1.0 / flux)
            flux = _
        elif flux.ndim == 1 and flux.shape[0] != 175:
            raise ValueError("Group structure must match EAF.")
        elif flux.ndim > 1:
            raise ValueError("The flux vector must be 0- or 1-dimensional.")
        if not np.all(flux >= 0.0):
            raise ValueError("Flux entries must be non-negative.")
        for ds in self.xs_cache.data_sources:
            ds.src_phi_g = flux
        self.xs_cache['phi_g'] = flux.sum()
        self._phi = flux

    def transmute(self, x, t=None, phi=None, log=None, tol=None):
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
                if part_nuc not in y_atoms:
                    y_atoms[part_nuc] = 0.0
                y_atoms[part_nuc] += part_adens * adens
        y = from_atom_frac(y_atoms)
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
            nuclide ids and values are float number densities for the coupled.

        """
        dest = self._get_destruction(nuc)
        A = np.zeros((1,1), float)
        A[0, 0] = -dest
        rootval = np.exp(-dest * self.t)
        partial = {nuc: rootval}
        self._traversal(nuc, A, partial)
        return partial

    def _get_destruction(nuc, decay=True):
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
        xs_chache = self.xs_cache
        sig_a = sigma_a(nuc, xs_cache=xs_cache)
        d = utils.from_barns(sig_a[0], 'cm2') * xs_cache['phi_g'][0]
        if decay:
            d += data.decay_const(nuc) 
        return d

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
            nuclide. Keys are nuclide names in integer (zzaaam) form. Values are
            number densities for the coupled nuclide in float format.  This is 
            modified in place.
        depth : int
            Current depth of traversal (root at 0). Should never be provided by user.

        """
        phi = self.phi
        temp = self.temp
        xs_cache = self.xs_cache
        if self.log is not None:
            self._log_tree(depth, nuc, 1.0)
        prod = {}
        # decay info
        lam = data.decay_const(nuc)
        decay_branches = {} if lam == 0 else self._decay_branches(nuc)
        for decay_child, branch_ratio in decay_branches.items():
            prod[decay_child] = lam * branch_ratio
        # reaction daughters
        daugh_dict = self._get_daughters(nuc)
        for rx in self.rxs:
            try:
                child = rxname.child(nuc, rx)
            except RuntimeError:
                continue
            child_xs = xs_cache[nuc, rx, temp][0]
            rr = utils.from_barns(child_xs, 'cm2') * phi  # reaction rate
            prod[child] = prod[child] + rr if child in prod else rr
        # Cycle production dictionary
        for child in prod_dict.keys():
            # Grow matrix
            B = _grow_matrix(A, prod_dict[child], _get_destruction(child, phi, table))
            # Create initial density vector
            n = B.shape[0]
            N0 = np.zeros((n,1))
            N0[0] = 1
            # Compute matrix exponential and dot with density vector
            eB = _matrix_exp(B, t)
            N_final = np.dot(eB, N0)
            # Log child
            if tree:
                _tree_log(depth+1, child, N_final[-1], tree)
            # Check against tolerance
            if _check_tol(N_final[-1], tol):
                # Continue traversal
                if tree is not None:
                    out = _traversal(child,B,phi,t,table,out,tol,tree,depth+1)
                else:
                    out = _traversal(child,B,phi,t,table,out,tol,tree,None)
            # On recursion exit or truncation, write data from this nuclide
            if child in out.keys():
                out[child] += N_final[-1]
            else:
                out[child] = N_final[-1]

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
        space = '   |'
        entry = "{spacing}--> {name} {N}\n".format(spacing=depth*space, 
                                                   numdens=numdens,
                                                   name=nucname.name(nuc))
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

#####################

def transmute_spatial(space, t_sim, tree = None, tol = 1e-7):
    """Transmutes a material into its daughters.

    Parameters
    ----------
    space : dictionary
        Input dictionary for the transmutation simulation.
        Keys are integers representing the volume ID.
        Values are tuples.
            The first entry in the tuple is a NumPy 1-dimensional array
                of floats representing the flux at the specified point.
            The second entry in the tuple is a the standard input
                dictionary used by the 'transmute' routine.
                (i.e. keys are zzaaam nuclides, values are initial
                    densities)
    t_sim : float
        Time to decay for.
    tree : File
        The file where the tree log should be written.
        tree should be None if a tree log is not desired.
    tol : float
        Tolerance level for chain truncation.
        Default tolerance level is 1e-7 for a root of unit density.

    Returns
    -------
    space_out : dictionary
        A dictionary containing the output from a multi-point simulation.
        Keys are integers representing the volume ID.
        Values are 'out' dictionaries (described below).
            out : dictionary
                A dictionary containing number densities for each
                nuclide after the simulation is carried out. Keys are
                nuclide names in integer (zzaaam) form. Values are
                number densities for the coupled nuclide in float format.
    """
    space_out = {}
    for volume in space.keys():
        phi = space[volume][0]
        inp = space[volume][1]
        out = transmute(inp, t_sim, phi, tree, tol)
        space_out[volume] = out
    return space_out




def _matrix_exp(A, t):
    """Takes the matrix exponential of the product At.

    Parameters
    ----------
    A : NumPy 2-dimensional array
        Coupled equation matrix.
    t : float
        Time required for transmutation simulation.

    Returns
    -------
    eA : NumPy 2-dimensional array
        Result after calculating the matrix exponential of At
    """
    eA = linalg.expm(A * t)
    return eA



def _grow_matrix(A, prod, dest):
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
    # Add row and column to current matrix
    B = np.append(A, np.zeros((1,n)), 0)
    B = np.append(B, np.zeros((n+1,1)), 1)
    # Update new matrix with provided data
    B[n,n-1] = prod
    B[n,n] = -dest
    return B


def _check_tol(N, tol):
    """Method to check if the current nuclide concentration exceeds the
    specified tolerance.

    Parameters
    ----------
    N : float
        The calculated final nuclide number density.
    tol : float
        The specified tolerance for the simulation.

    Returns
    -------
    fail : Boolean
        False if the final nuclide density is less than the tolerance.
        True if the final nuclide density is not less than the tolerance.
    """
    fail = N > tol
    return fail



