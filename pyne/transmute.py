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


def decay(nuc, phi, t_sim, tol):
    """Decays a material into its daughters.

    Parameters
    ----------
    nuc : nucname
        Name of the nuclide in decay.
    phi : NumPy array of floats
        Neutron flux vector of length 175.
    t_sim : float
        Time to decay for.
    tol : float
        Tolerance level for chain truncation.

    Returns
    -------
    decay_nuc : NumPy array of nucnames
        The daughters.
    """
    # Check parameters
    # Convert nuc to zzaaam
    nuc = nucname.zzaaam(nuc)
    # Check length of phi
    if not(phi.size == 175):
        sys.exit('Incorrect phi dimension for FENDL data.')
    A = _create_decay_matrix(nuc)
    eA = _solve_decay_matrix(A)
    data_table.close()


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
    eA = sp.linalg.expm(A * t)
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
        daugh_dict[daugh] = xs
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
    xs_total = np.zeros((175))
    for nuc in rxn_dict.keys():
        xs_total = xs_total + rxn_dict[nuc]
    decay_const = data.decay_const(nuc)
    d = decay_const + sum(xs_total*phi)
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


def _check_tol(N, tol)
    """Method to check if the current nuclide concentration exceeds the
    specified tolerance.

    Parameters
    ----------
    N : NumPy 1-dimensional array
        The calculated vector of nuclide number densities
    tol : float
        The specified tolerance for the simulation.

    Returns
    -------
    fail : Boolean
        False if the final nuclide density is less than the tolerance.
        True if the final nuclide density is not less than the tolerance.
    """
    fail = N_final[-1] > tol
    return fail


def _traversal(nuc, A, phi, t, N_ini, out, tol):
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
    
    Returns
    -------
    out : dictionary
        A dictionary containing the final recorded number densities for each
        nuclide. Keys are nuclide names in integer (zzaaam) form. Values are
        number densities for the coupled nuclide in float format.
    """
    # Lookup decay constant of current nuclide
    lam = data.decay_const(nuc)
    # Lookup decay products and reaction daughters
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
            prod_dict[decay_daugh] += sum(phi * daugh_dict[decay_daugh])
        else:
            prod_dict[decay_daugh] = sum(phi * daugh_dict[decay_daugh])
    # Cycle production dictionary
    for child in prod_dict.keys():
        # Create initial density vector
        shape = B.shape
        n = shape[0]
        N0 = np.zeros((n,1))
        N0[0] = N_ini
        # Grow matrix
        B = _grow_matrix(A, prod_dict[child], _get_destruction(child, phi))
        # Compute matrix exponential and dot with density vector
        eB = _matrix_exp(B, t)
        N_final = np.dot(eB, N0)
        # Check against tolerance
        if _check_tol(N_final, tol):
            # Continue traversal
            out = _traversal(child, B, phi, t, N_ini, out, tol)
        # On recursion exit or truncation, write data from this nuclide
        if child in out.keys():
            out[child] += N_final[-1]
        else:
            out[child] = N_final[-1]
    # Return final output dictionary
    return out

