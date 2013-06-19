import numpy as np
import scipy as sp
import tables as tb
import math
import sys

from pyne import data
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material
from pyne.xs.data_source import EAF_RX
from scipy import linalg
from tables import *


def transmute_core(nuc, t_sim, phi, tree = None, tol = 1e-7):
    """Core method to transmute a material into its daughters.
    This method assumes that the initial nuclide has unit density.

    Parameters
    ----------
    nuc : nucname
        Integer representation of nuclide to be transmuted.
    t_sim : float
        Time to decay for.
    phi : NumPy 1-dimensional array of floats
        Neutron flux vector.
        If phi is None, the flux vector is set to zero.
    tree : File
        The file where the tree log should be written.
        tree should be None if a tree log is not desired.
    tol : float
        Tolerance level for chain truncation.
        Default tolerance level is 1e-7 for a root of unit density.

    Returns
    -------
    out : dictionary
        A dictionary containing number densities for each nuclide after
        the simulation is carried out. Keys are nuclide names in integer
        (zzaaam) form. Values are number densities for the coupled
        nuclide in float format.
    """
    out = {}
    phi = _check_phi(phi)
    dest = _get_destruction(nuc, phi)
    A = np.zeros((1,1))
    A[0,0] = -dest
    rootVal = math.exp(-dest * t_sim)
    out = {nuc : rootVal}
    out = _traversal(nuc, A, phi, t_sim, out, tol, tree, depth = None)
    return out


def transmute(inp, t_sim, phi, tree = None, tol = 1e-7):
    """Transmutes a material into its daughters.

    Parameters
    ----------
    inp : dictionary
        Input dictionary for the transmutation simulation.
        Keys are nuclides in integer (zzaaam) format.
        Values are corresponding number densities represented by floats.
    t_sim : float
        Time to decay for.
    phi : NumPy 1-dimensional array of floats
        Neutron flux vector.
        If phi is None, the flux vector is set to zero.
    tree : File
        The file where the tree log should be written.
        tree should be None if a tree log is not desired.
    tol : float
        Tolerance level for chain truncation.
        Default tolerance level is 1e-7 for a root of unit density.

    Returns
    -------
    out : dictionary
        A dictionary containing number densities for each nuclide after
        the simulation is carried out. Keys are nuclide names in integer
        (zzaaam) form. Values are number densities for the coupled
        nuclide in float format.
    """
    # Properly format phi
    phi = _check_phi(phi)
    out = {}
    for nuc in inp.keys():
        # Find output for root of unit density
        out_partial = transmute_core(nuc, t_sim, phi, tree, tol)
        # Scale all output by actual nuclide density
        for part in out_partial.keys():
            out_partial[part] = out_partial[part] * inp[nuc]
    return out


def transmute_spatial(space, t_sim, tree = None, tol = 1e-7):
    """Transmutes a material into its daughters.

    Parameters
    ----------
    space : dictionary
        Input dictionary for the transmutation simulation.
        Keys are float triples corresponding to the xyz location in space.
        Values are tuples.
            The first entry in the tuple is a NumPy 1-dimensional array
                of floats representing the flux at the specified point.
            The second entry in the tuple is a the standard input
                dictionary used by the 'transmute' routine.
                (i.e. keys are zzaaam nuclides, values are initial
                    densities)
    t_sim : float
        Time to decay for.
    phi : NumPy 1-dimensional array of floats
        Neutron flux vector.
        If phi is None, the flux vector is set to zero.
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
        Keys are float triples corresponding to the xyz location in space.
        Values are 'out' dictionaries (described below).
            out : dictionary
                A dictionary containing number densities for each
                nuclide after the simulation is carried out. Keys are
                nuclide names in integer (zzaaam) form. Values are
                number densities for the coupled nuclide in float format.
    """
    space_out = {}
    for point in space.keys():
        phi = space[point][0]
        inp = space[point][1]
        out = transmute(inp, t_sim, phi, tree, tol)
        space_out[point] = out
    return space_out


class Nuclide(IsDescription):
    """Class to describe columns of hdf5 output table"""
    size = 16
    name = StringCol(size)
    zzaaam = StringCol(size)
    density = FloatCol()


def write_hdf5(h5file, parentGroup, out, title = None):
    """Method to write contents of an output dictionary generated by
    transmute() to a specified group of an hdf5 file.

    Parameters
    ----------
    h5file : PyTables file handle
        The PyTables file representation of the hdf5 file that should
        be written to.
    parentGroup : String
        The String representation of the parent group that the output
        should be written under.
    out : dictionary
        A dictionary containing number densities for each nuclide after
        the simulation is carried out. Keys are nuclide names in integer
        (zzaaam) form. Values are number densities for the coupled
        nuclide in float format.
    title : String
        Optionally, a title to assign to the table being written

    Returns
    -------
    None
        This method writes to a file and does not return any
        information.
    """
    if title is None:
        table = h5file.createTable(parentGroup,'transmutation',Nuclide)
    else:
        table = h5file.createTable(parentGroup, title, Nuclide)
    nuc = table.row
    for key in out.keys():
        nuc['name'] = nucname.name(key)
        nuc['zzaaam'] = key
        nuc['density'] = out[key]
        nuc.append()
    h5file.flush()
    return None


def write_space_hdf5(h5file, parentGroup, space_out):
    """Method to write contents of an output dictionary generated by
    transmute() to a specified group of an hdf5 file.

    Parameters
    ----------
    h5file : PyTables file handle
        The PyTables file representation of the hdf5 file that should
        be written to.
    parentGroup : String
        The String representation of the parent group that the output
        should be written under.
    space_out : dictionary
        A dictionary containing the output from a multi-point simulation.
        Keys are float triples corresponding to the xyz location in space.
        Values are 'out' dictionaries (described below).
            out : dictionary
                A dictionary containing number densities for each
                nuclide after the simulation is carried out. Keys are
                nuclide names in integer (zzaaam) form. Values are
                number densities for the coupled nuclide in float format.

    Returns
    -------
    None
        This method writes to a file and does not return any
        information.
    """
    node = h5file.createGroup(parentGroup,'transmutation', \
                                'Multi-point transmutation')
    for point in space_out.keys():
        out = space_out[point]
        write_hdf5(h5file, node, out, 'Point: ' + str(point))
    h5file.flush()
    return None
            

def _check_phi(phi):
    """Ensures that the flux vector phi is correctly formatted.

    Parameters
    ----------
    phi : NumPy 1-dimensional array
        Phi may be either correctly formatted or None.
        When the value of phi is None, a vector of zero flux is used.

    Returns
    -------
    phi : NumPy 1-dimensional array
        Phi will be returned with a shape of (numEntries,1).
    """
    numEntries = 175
    if phi is None:
        phi = np.zeros((numEntries,1))
        return phi
    n = phi.shape[0]
    if phi.ndim != 2:
        # Throw exception for incorrect shape
        pass
    if n != numEntries:
        # Throw exception for incorrect number of entries
        pass
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
    eA = linalg.expm2(A * t)
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
        NOTE
            Cross sections have been converted from units of b to units
            of cm^2.
    """
    numEntries = 175
    barn_cm2 = 1e-24
    fissionMT = 180
    daugh_dict = {}
    # Remove fission MT# (cannot handle)
    if fissionMT in EAF_RX:
        EAF_RX.remove(fissionMT)
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
        # Convert from barns to cm^2
        xs = all_xs[i] * barn_cm2
        daugh_dict[daugh] = xs.reshape((numEntries,1))
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


def _get_destruction(nuc, phi, addDecay = True):
    """Returns the destruction rate of the nuclide.

    Parameters
    ----------
    nuc : nucname
        Name of the nuclide in question
    phi : NumPy 1-dimensional array
        Flux vector for use in simulation. The vector should be 175 entries
        in length for proper usage with EAF data.
    addDecay : boolean
        True if the decay constant should be added to the returned value.
        False if only destruction from neutron reactions should be
            considered.

    Returns
    -------
    d : float
        Destruction rate of the nuclide.
    """
    numEntries = 175
    nuc = nucname.zzaaam(nuc)
    rxn_dict = _get_daughters(nuc)
    xs_total = np.zeros((numEntries,1))
    for key in rxn_dict.keys():
        xs_total += rxn_dict[key]
    if addDecay:
        d = data.decay_const(nuc) + np.sum(xs_total*phi)
    else:
        d = np.sum(xs_total*phi)
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
    spacing = depth * '   |'
    name = nucname.name(nuc)
    Nstr = str(N)
    entry = spacing + '--> ' + name + ' (' + Nstr + ')\n'
    tree.write(entry)
    return None


def _traversal(nuc, A, phi, t, out, tol, tree, depth = None):
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
            out = _traversal(child, B, phi, t, out, tol, tree, depth+1)
        # On recursion exit or truncation, write data from this nuclide
        if child in out.keys():
            out[child] += N_final[-1]
        else:
            out[child] = N_final[-1]
    # Return final output dictionary
    return out

