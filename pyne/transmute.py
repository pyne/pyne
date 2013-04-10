import numpy as np
import scipy as sp
import tables as tb
import sys
import os

from pyne import data
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material
from pyne.xs.data_source import EAF_RX
from scipy import linalg


def transmute(inp, t_sim, tol, tree, phi = None, filename = None):
    """Transmutes a material into its daughters.

    Parameters
    ----------
    inp : dictionary
        Input dictionary for the transmutation simulation.
        Keys are nuclides in integer (zzaaam) format.
        Values are corresponding number densities represented by floats.
    t_sim : float
        Time to decay for.
    tol : float
        Tolerance level for chain truncation.
    tree : Boolean
        True if a tree output file is desired.
        False if a tree output file is not desired.
    phi : NumPy 1-dimensional array of floats
        Neutron flux vector.
        If phi is None, the flux vector is set to zero.
        If phi is less than 175 entries in length, zeros will be added
        until it contains 175 entries.
    filename : String
        Name of file to write tree log to.
        Must be provided if 'tree' is True.

    Returns
    -------
    out : dictionary
        A dictionary containing number densities for each nuclide after
        the simulation is carried out. Keys are nuclide names in integer
        (zzaaam) form. Values are number densities for the coupled
        nuclide in float format.
    """
    # Check length of phi
    if phi is None:
        phi = np.zeros((175, 1))
    else:
        phi = _format_phi(phi)
    out = {}
    for nuc in inp.keys():
        A = np.zeros((1,1))
        dest = _get_destruction(nuc,phi)
        A[0,0] = -dest
        N_ini = inp[nuc]
        out = _traversal(nuc, A, phi, t_sim, N_ini, out, tol, tree, filename,\
                            True)
    return out


def _format_phi(phi):
    """Ensures that the flux vector phi is correctly formatted.

    Parameters
    ----------
    phi : NumPy 1-dimensional array
        Phi may be of various acceptable formats.

    Returns
    -------
    phi : NumPy 1-dimensional array
        Phi will be returned with a shape of (175,1).
    """
    n = phi.shape[0]
    if phi.ndim == 1:
        phi = phi.reshape((n,1))
    if n < 175:
        rem = 175 - n
        app = np.zeros((rem,1))
        phi = np.append(phi, app, 0)
    return phi


def _create_decay_matrix(nucs):
    """Creates a NumPy matrix representing the nuclide list nucs.

    Parameters
    ----------
    nucs : NumPy array of nucnames
        Names of nuclides to build decay matrix around.

    Returns
    -------
    A : NumPy matrix of floats
        The decay matrix.
    """
    nnucs = len(nucs)
    nucsrange = np.arange(nnucs)
    A = np.zeros((nnucs, nnucs), dtype=float)
    A[nucsrange, nucsrange] = [-data.decay_const(nuc) for nuc in nucs]
    return A


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


def _get_daughters(nuc):
    """Returns a dictionary that contains the neutron-reaction daughters of
    nuc as keys to the 175-group neutron cross sections for that daughter's
    reaction.

    Parameters
    ----------
    nuc : nucname
        Name of parent nuclide to get daughters of.

    Returns
    -------
    daugh_dict : dictionary
        Keys are the neutron-reaction daughters of nuc in zzaaam format.
        Values are a NumPy array containing the EAF cross section data.
            (all Values should have size 175)
    """
    daugh_dict = {}
    # Remove fission MT# (cannot handle)
    if '180' in EAF_RX:
        EAF_RX.remove('180')
    # Open nuc_data.h5
    with tb.openFile(nuc_data, 'r') as f:
        # Set working node that contains EAF cross sections
        node = f.root.neutron.eaf_xs.eaf_xs
        cond = "(nuc_zz == {0})".format(nuc)
        daughters = [row['daughter'] for row in node.where(cond)]
        all_xs = [np.array(row['xs']) for row in node.where(cond)]
        all_rx = [row['rxnum'] for row in node.where(cond)]
    for i in range(len(daughters)):
        if all_rx[i] not in EAF_RX:
            continue
        daugh = _convert_eaf(daughters[i])
        xs = all_xs[i]
        daugh_dict[daugh] = xs.reshape((175,1))
    return daugh_dict


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
    # Remove possible space from string
    daugh = daugh.replace(' ', '')
    # Check for 'G' suffix
    daugh = daugh.replace('G', '')
    # Check for metastable suffix
    if daugh.endswith('M1') or daugh.endswith('M2'):
        # Convert appropriately
        parts = daugh.rsplit('M',1)
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


def _get_destruction(nuc, phi):
    """Returns the destruction rate of the nuclide.

    Parameters
    ----------
    nuc : nucname
        Name of the nuclide in question
    phi : NumPy 1-dimensional array
        Flux vector for use in simulation. The vector should be 175 entries
        in length for proper usage with EAF data.

    Returns
    -------
    d : float
        Destruction rate of the nuclide.
    """
    nuc = nucname.zzaaam(nuc)
    rxn_dict = _get_daughters(nuc)
    xs_total = np.zeros((175,1))
    for key in rxn_dict.keys():
        xs_total += rxn_dict[key]
    d = data.decay_const(nuc) + np.sum(xs_total*phi)
    return d


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


def _tree_log(depth, nuc, N, filename, new = None):
    """Logging method to track path of _traversal.

    Parameters
    ----------
    depth : integer
        Current depth of traversal (root at 0).
    nuc : nucname
        Current nuclide in traversal.
    N : float
        Current density of nuc.
    filename : String
        Name of file to write tree log to.
    new : boolean
        True if a new file should be created or existing file should be
            overwritten.
        False or None if the current filename should be appended to.

    Returns
    -------
    None
    """
    spacing = depth * '----'
    name = nucname.name(nuc)
    Nstr = str(N)
    entry = spacing + '> ' + name + ' (' + Nstr + ')\n'
    if new:
        with open(filename, 'w') as f:
            f.write(entry)
    else:
        with open(filename, 'a') as f:
            f.write(entry)
    return None


def _traversal(nuc, A, phi, t, N_ini, out, tol, tree, filename = None, \
                first = None, depth = None):
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
    N_ini : float
        Number density of initial nuclide at root of transmutation tree.
    out : dictionary
        A dictionary containing the final recorded number densities for each
        nuclide. Keys are nuclide names in integer (zzaaam) form. Values are
        number densities for the coupled nuclide in float format.
    tol : float
        Tolerance level to reference for tree truncation.
    tree : Boolean
        True if a tree output file is desired.
        False if a tree output file is not desired.
    filename : String
        Name of file to write tree log to.
        Must be provided if 'tree' is True.
    first : boolean
        True if this is the first traversal for the problem,
            False or None otherwise.
        This argument is only required if a tree log is desired.
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
    if depth is None and tree:
        if filename is None:
            # Filename not provided, throw exception
            pass
        depth = 0
        if first:
            _tree_log(depth, nuc, N_ini, filename, True)
        else:
            _tree_log(depth, nuc, N_ini, filename)
    # Lookup decay constant of current nuclide
    lam = data.decay_const(nuc)
    # Lookup decay products and reaction daughters
    if lam == 0:
        decay_dict = {}
    else:
        decay_dict = _get_decay(nuc)
    daugh_dict = _get_daughters(nuc)
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
        B = _grow_matrix(A, prod_dict[child], _get_destruction(child, phi))
        # Create initial density vector
        n = B.shape[0]
        N0 = np.zeros((n,1))
        N0[0] = N_ini
        # Compute matrix exponential and dot with density vector
        eB = _matrix_exp(B, t)
        N_final = np.dot(eB, N0)
        # Log child
        if tree:
            _tree_log(depth+1, child, N_final[-1], filename)
        # Check against tolerance
        if _check_tol(N_final[-1], tol):
            # Continue traversal
            out = _traversal(child, B, phi, t, N_ini, out, tol, tree,\
                                filename, False, depth+1)
        # On recursion exit or truncation, write data from this nuclide
        if child in out.keys():
            out[child] += N_final[-1]
        else:
            out[child] = N_final[-1]
    # Return final output dictionary
    return out

