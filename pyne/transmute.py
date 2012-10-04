import numpy as np
import scipy as sp

from pyne.material import Material
from pyne import data


def decay(mat, t):
    """Decays a material into its daughters.

    Parameters
    ----------
    mat: Material 
        Material to decay.
    t : float
        Time to decay for.

    Returns
    -------
    decay_mat : Material
        The daughters.
    """
    nucvec = mat.keys()
    A = _create_decay_matrix(nucvec)
    eA = _solve_decay_matrix(A)


def _create_decay_matrix(nucs):
    nnucs = len(nucs)
    nucsrange = np.arange(nnucs)
    A = np.zeros((nnucs, nnucs), dtype=float)
    A[nucsrange, nucsrange] = [-data.decay_const(nuc) for nuc in nucs]
    return A


def _solve_decay_matrix(A):
    eA = sp.linalg.expm(A)
    return eA
