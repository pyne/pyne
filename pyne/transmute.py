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


def _solve_decay_matrix(A):
    eA = sp.linalg.expm(A)
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

def _get_destruction(nuc, phi, rxn_dict = None):
    """Returns the destruction rate of the nuclide.

    Parameters
    ----------
    nuc : nucname
        Name of the nuclide in question

    rxn_dict : dictionary
        Contains the dictionary of neutron reactions for nuc that follows the
        structure of _get_daughters() output. This variable avoids a second
        lookup of the information if the user has already called 
        _get_daughters().

    Returns
    -------
    d : float
        Destruction rate of the nuclide.
    """
    nuc = nucname.zzaaam(nuc)
    if rxn_dict is None:
        rxn_dict = _get_daughters(nuc)
    xs_total = np.zeros((175))
    for nuc in rxn_dict.keys():
        xs_total = xs_total + rxn_dict[nuc]
    decay_const = data.decay_const(nuc)
    d = decay_const + sum(xs_total*phi)
    return d


