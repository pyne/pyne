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
    # Store current matrix size
    shape = A.shape
    n = shape[0]
    # Find decay constant of current nuclide
    lam = data.decay_const(nuc)
    # Cycle decay children
    decay_dict = _get_decay(nuc)
    for decay_child in decay_dict.keys():
        # Lookup branch ratio
        branch_rat = decay_dict[key]
        # Find destruction rate of decay_child for appending to matrix
        dest = _get_destruction(decay_child, phi)
        # Add row and column to matrix
        B = A
        B = np.append(A, np.zeros((1,n)), 0)
        B = np.append(B, np.zeros((n+1,1)), 1)
        # Create initial density vector
        N0 = np.zeros((n+1,1))
        N0[0] = N_ini
        # Determine production rate of decay_child
        prod = lam * branch_rat
        # Update new matrix B
        B[n,n-1] = prod
        B[n,n] = -dest
        # Calculate density of decay_child
        eB = _matrix_exp(B, t)
        N_final = np.dot(eB, N0)
        # Check against tolerance
        if N_final[-1] > tol:
            # Continue traversal
            out = _traversal(decay_child, B, phi, t, N_ini, out, tol)
        # On recursion exit or truncation, write data from this nuclide
        if decay_child in out.keys():
            # If already in output dictionary, increment instead
            out[decay_child] = out[decay_child] + N_final[-1]
        else:
            out[decay_child] = N_final[-1]
       
    # Cycle neutron reaction daughters
    daugh_dict = _get_daughters(nuc)
    for daugh in daugh_dict.keys():
        # Lookup appropriate cross section for daugh
        xs = daugh_dict[daugh]
        # Find destruction rate of daugh
        dest = _get_destruction(daugh, phi)
        # Add row and column to matrix
        B = A
        B = np.append(A, np.zeros((1,n)), 0)
        B = np.append(B, np.zeros((n+1,1)), 1)
        # Create initial density vector
        N0 = np.zeros((n+1,1))
        N0[0] = N_ini
        # Determine production rate of daugh
        prod = sum(xs*phi)
        # Update new matrix B
        B[n,n-1] = prod
        B[n,n] = dest
        # Calculate density of daugh
        eB = _matrix_exp(B, t)
        N_final = np.dot(eB,N0)
        # Check against tolerance
        if N_final[-1] > tol:
            # Continue traversal
            out = _traversal(daugh, B, phi, t, N_ini, out, tol)
        # On recursion exit or truncation, write data from this nuclide
        if daugh in out.keys():
            # If already in output dictionary, increment density
            out[daugh] = out[daugh] + N_final[-1]
        else:
            out[daugh] = N_final[-1]

    return out

