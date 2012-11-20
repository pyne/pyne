import numpy as np
import scipy as sp
import tables as tb
import sys

from pyne import data
from pyne import nucname
from pyne import nuc_data
from pyne.material import Material
from pyne.xs.cache import xs_cache
from pyne.dbgen import eaf#, BASIC_FILTERS


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
    # Convert nuc to zzaaam
    nuc = nucname.zzaaam(nuc)
    # Check length of phi
    if not(phi.size == 175):
        sys.exit("Incorrect phi dimension for FENDL data.")
    A = _create_decay_matrix(nuc)
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


def _import_eaf_data():
    """Creates hdf5 table "nuc_data" in local directory from FENDL data.

    Location of EAF data hard-coded for CNERG machines.
    
    Method may be deprecated; hdf5 table already exists, though location has
    been hard-coded again.

    Returns
    -------
    eaf_table : PyTables table
        hdf5 table format of FENDL data.
    """

    eaf.make_eaf_table("nuc_data", \
                   "/filespace/groups/cnerg/opt/FENDL2.0-A/fendlg-2.0_175")


def _get_daughters(nuc):
    """Returns an array of nuclides which are daughters of nuc in a numpy
    array.

    Returns
    -------
    daughters : list
        Daughter nuclides of nuc in human-readable format.
    """

    eaf_table = tb.openFile('../build_nuc_data/prebuilt_nuc_data.h5')
    daughters = [row['daughter'] for row in \
        eaf_table.root.neutron.eaf_xs.eaf_xs.where('nuc_zz == nuc')]
#Problem with daughter formatting (e.g. '.. 0')
    #for i in range(len(daughters)):
        #daughters[i] = nucname.zzaaam(daughters[i].replace(' ',''))
    eaf_table.close()
    return daughters
