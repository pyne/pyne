"""This module implements an ALARA-like chain-based transmutation solver.
"""

import numpy as np
from scipy import linalg

from pyne import data
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material, from_atom_frac
from pyne.xs.data_source import EAF_RX, NullDataSource, EAFDataSource
from pyne.xs.cache import XSCache
from pyne.xs.channels import sigma_a

class Transmuter(object):
    """A class for transmuting materials using an ALARA-like chain solver."""

    def __init__(self, t=0.0, phi=0.0, tol=1e-7, log=None):
        """Parameters
        ----------
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
        """
        self.t = t
        self._phi = None
        self.phi = phi
        self.log = log
        self.tol = tol
        
        eafgs = EAFDataSource().src_group_stuct
        gs = np.array([eafgs[0], eafgs[-1]])
        self.xs_cache = XSCache(group_struct=gs, 
                                data_source_classes=(EAFDataSource, NullDataSource))

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, flux):
        """Ensures that the flux is correctly formatted."""
        flux = np.asarray(flux)
        if flux.ndim > 1:
            raise ValueError("The flux vector must be 0- or 1-dimensional.")
        if flux.ndim == 1 and flux.shape[0] != 175:
            raise ValueError("Group structure must match EAF.")
        if not np.all(flux >= 0.0):
            raise ValueError("Flux entries must be non-negative.")
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
        partial = {}
        dest = self._get_destruction(nuc)
        A = np.zeros((1,1))
        A[0,0] = -dest
        rootval = np.exp(-dest * self.t)
        partial = {nuc: rootval}
        partial = _traversal(nuc, A, phi, t_sim, table, out, tol, tree, depth = None)
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
        rxn_dict = _get_daughters(nuc)
        xs_total = np.zeros((eaf_numEntries, 1))
        for key in rxn_dict.keys():
            xs_total += rxn_dict[key]
        d = np.sum(xs_total * phi)
        if decay:
            d += data.decay_const(nuc) 
        return d

    def _get_daughters(self, nuc):
        """Returns a dictionary that contains the neutron-reaction daughters of
        nuc as keys to the 175-group neutron cross sections for that daughter's
        reaction.

        Parameters
        ----------
        nuc : int
            ID of mother nuclide to of which to find daughters.

        Returns
        -------
        daugh_dict : dictionary
            Keys are the neutron-reaction daughters of nuc in zzaaam format.
            Values are a NumPy array containing the EAF cross section data.
            (all Values should have size 175)
            NOTE
                Cross sections have been converted from units of b to units
                of cm^2.
        """
        eaf_numEntries = 175
        barn_cm2 = 1e-24
        fissionMT = 180
        daugh_dict = {}
        # Remove fission MT# (cannot handle)
        if fissionMT in EAF_RX:
            EAF_RX.remove(fissionMT)
        # Set working node that contains EAF cross sections
        node = table.root.neutron.eaf_xs.eaf_xs
        cond = "(nuc_zz == {0})".format(nuc)
        daughters = [row['daughter'] for row in node.where(cond)]
        all_xs = [np.array(row['xs']) for row in node.where(cond)]
        all_rx = [row['rxnum'] for row in node.where(cond)]
        for i in range(len(daughters)):
            if all_rx[i] not in EAF_RX:
                continue
            daugh = _convert_eaf(daughters[i])
            # Convert from barns to cm^2
            xs = all_xs[i] * barn_cm2
            daugh_dict[daugh] = xs.reshape((eaf_numEntries, 1))
        return daugh_dict




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



def _convert_eaf(daugh):
    """Returns the zzaaam format of a daugh string in parsed EAF format.

    Parameters
    ----------
    daugh : String
        String representation of a daughter in EAF format.

    Returns
    -------
    daugh_conv : nucname
        Name of daugh in zzaaam format appropriate for PyNE.
    """
    singleSpace = ' '
    noSpace = ''
    letterG = 'G'
    meta1 = 'M1'
    meta2 = 'M2'
    letterM = 'M'
    # Remove possible space from string
    daugh = daugh.replace(singleSpace, noSpace)
    # Check for 'G' suffix
    if daugh.endswith(letterG):
        daugh_noG  = daugh.rsplit(letterG, 1)[0]
        daugh_conv = nucname.zzaaam(daugh_noG)
    # Check for metastable suffix
    elif daugh.endswith(meta1) or daugh.endswith(meta2):
        # Convert appropriately
        parts = daugh.rsplit(letterM ,1)
        daugh_conv = nucname.zzaaam(parts[0]) + int(parts[1])
    else:
        daugh_conv = nucname.zzaaam(daugh)
    return daugh_conv


def _get_decay(nuc):
    """Returns a dictionary that contains the decay children of nuc as keys
    to the branch ratio of that child's decay process.

    Parameters
    ----------
    nuc : nucname
        Name of parent nuclide to get decay children of.

    Returns
    -------
    decay_dict : dictionary
        Keys are decay children of nuc in zzaaam format.
        Values are the branch ratio of the decay child.
    """
    decay_dict = {}
    nuc = nucname.zzaaam(nuc)
    children = data.decay_children(nuc)
    for child in children:
        branch = data.branch_ratio(nuc,child)
        decay_dict[child] = branch
    return decay_dict



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


def _tree_log(depth, nuc, N, tree):
    """Logging method to track path of _traversal.

    Parameters
    ----------
    depth : integer
        Current depth of traversal (root at 0).
    nuc : nucname
        Current nuclide in traversal.
    N : float
        Current density of nuc.
    tree : File
        File to write tree log to.

    Returns
    -------
    None
        This method only writes to the File "tree".
    """
    # Don't print a zero density.
    if N == 0:
        return None
    space = '   |'
    arrow = '--> '
    blank = ' '
    newline = '\n'
    spacing = depth * space
    name = nucname.name(nuc)
    Nstr = str(N)
    entry = spacing + arrow + name + blank + Nstr + newline
    tree.write(entry)
    tree.flush()
    return None


def _traversal(nuc, A, phi, t, table, out, tol, tree, depth = None):
    """Nuclide transmutation traversal method.

    This method will traverse the reaction tree recursively, using a DFS
    algorithm. On termination, the method will return all number densities
    after a given time that are a result of the starting nuclide.

    Parameters
    ----------
    nuc : nucname
        Name of the active nuclide for the traversal.
    A : NumPy 2-dimensional array
        Current state of the coupled equation matrix.
    phi : NumPy 1-dimensional array
        Flux vector for use in simulation. The vector should be 175 entries
        in length for proper usage with EAF data.
    t : float
        Time at which to evaluate transmutation events.
    table : hdf5 file handle
        The nuc_data hdf5 table handle.
    out : dictionary
        A dictionary containing the final recorded number densities for each
        nuclide. Keys are nuclide names in integer (zzaaam) form. Values are
        number densities for the coupled nuclide in float format.
    tol : float
        Tolerance level to reference for tree truncation.
    tree : Boolean
        True if a tree output file is desired.
        False if a tree output file is not desired.
    depth : integer
        Current depth of traversal (root at 0).
        Should never be provided by user.

    Returns
    -------
    out : dictionary
        A dictionary containing the final recorded number densities for each
        nuclide. Keys are nuclide names in integer (zzaaam) form. Values are
        number densities for the coupled nuclide in float format.
    """
    # Log initial nuclide
    if depth is None and tree is not None:
        depth = 0
        _tree_log(depth, nuc, 1, tree)
    # Lookup decay constant of current nuclide
    lam = data.decay_const(nuc)
    # Lookup decay products and reaction daughters
    if lam == 0:
        decay_dict = {}
    else:
        decay_dict = _get_decay(nuc)
    daugh_dict = _get_daughters(nuc, table)
    # Initialize production rate dictionary
    prod_dict = {}
    # Cycle decay children
    for decay_child in decay_dict.keys():
        prod_dict[decay_child] = lam * decay_dict[decay_child]
    # Cycle reaction daughters
    for decay_daugh in daugh_dict.keys():
        # Increment current production rate if already in dictionary
        if decay_daugh in prod_dict.keys():
            prod_dict[decay_daugh] += np.sum(phi * daugh_dict[decay_daugh])
        else:
            prod_dict[decay_daugh] = np.sum(phi * daugh_dict[decay_daugh])
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
    # Return final output dictionary
    return out

