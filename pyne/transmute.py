import numpy as np
import scipy as sp
import tables as tb

from pyne import data
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material
#from pyne.xs.cache import xs_cache
from pyne.dbgen import eaf#, BASIC_FILTERS


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

def import_eaf_data():
    """Creates numpy array by parsing EAF data.

    Location of EAF data hard-coded for CNERG machines

    Returns
    -------
    eaf_table : PyTables table

    """

    #eaf_array = eaf.parse_eaf_xs('/filespace/groups/cnerg/opt/FENDL2.0-A/fendlg-2.0_175')
    eaf.make_eaf_table("nuc_data", \
                   "/filespace/groups/cnerg/opt/FENDL2.0-A/fendlg-2.0_175")
    #return eaf_array

def _get_daughters(nuc):
    eaf_table = tb.openFile('nuc_data')
    daughters = [row['daughter'] for row in \
        eaf_table.root.neutron.eaf_xs.eaf_xs.where('nuc_zz == nuc')]
    return daughters
