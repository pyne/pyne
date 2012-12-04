import numpy as np
import scipy as sp
import tables as tb
import sys
import os

from pyne import data
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material
from pyne.dbgen import eaf


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
    # Fetch nuc_data.h5
    data_table = tb.openFile(nuc_data)
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


def _get_daughters(nuc, data_table):
    """Returns an array of nuclides which are daughters of nuc in a numpy
    array.

    Parameters
    ----------
    nuc : nucname
        Name of parent nuclide to get daughters of.
    data_table : PyTables table
        nuc_data.h5 file

    Returns
    -------
    daughters : list
        Daughter nuclides of nuc in human-readable format.
    """
    daughters = [row['daughter'] for row in \
        data_table.root.neutron.eaf_xs.eaf_xs.where('nuc_zz == nuc')]
# Problem with daughter formatting (e.g. '.. 0')
# This code would convert the names from human to zzaaam
    #for i in range(len(daughters)):
        #daughters[i] = nucname.zzaaam(daughters[i].replace(' ',''))
    return daughters



