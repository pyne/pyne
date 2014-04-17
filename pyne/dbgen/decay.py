"""This module provides a way to grab and store raw data for radioactive decay."""
from __future__ import print_function, division
import os
import glob
try:
    import urllib.request as urllib
except ImportError:
    import urllib
from zipfile import ZipFile

import numpy as np
import tables as tb

from pyne import ensdf
from pyne.dbgen.api import BASIC_FILTERS


def grab_ensdf_decay(build_dir=""):
    """Grabs the ENSDF decay data files
    if not already present.

    Parameters
    ----------
    build_dir : str
        Major directory to place html files in. 'ENSDF/' will be appended.
    """
    # Add ENSDF to build_dir
    build_dir = os.path.join(build_dir, 'ENSDF')
    try:
        os.makedirs(build_dir)
    except OSError:
        pass

    # Grab ENSDF files and unzip them.
    iaea_base_url = 'http://www-nds.iaea.org/ensdf_base_files/2013-October/'
    cf_base_url = 'http://data.pyne.io/'
    ensdf_zip = ['ensdf_131009_099.zip', 'ensdf_131009_199.zip', 'ensdf_131009_294.zip', ]

    for f in ensdf_zip:
        fpath = os.path.join(build_dir, f)
        if f not in os.listdir(build_dir):
            print("  grabbing {0} and placing it in {1}".format(f, fpath))
            urllib.urlretrieve(iaea_base_url + f, fpath)

            if os.path.getsize(fpath) < 1048576:
                print("  could not get {0} from IAEA; trying mirror".format(f))
                os.remove(fpath)
                urllib.urlretrieve(cf_base_url + f, fpath)

        # not using ZipFile context manager (with statement for Python 2.6)
        try:
            zf = ZipFile(fpath)
            for name in zf.namelist():
                if not os.path.exists(os.path.join(build_dir, name)):
                    print("    extracting {0} from {1}".format(name, f))
                    zf.extract(name, build_dir)
        finally:
            zf.close()


level_dtype = np.dtype([
    ('nuc_id', int),
    ('rx_id', np.uint32),
    ('half_life', float),
    ('level', float),
    ('branch_ratio', float),
    ('metastable', int),
    ('special', 'S1'),
])

decay_dtype = np.dtype([
    ('parent', int),
    ('child', int),
    ('decay', np.uint32),
    ('half_life', float),
    ('half_life_error', float),
    ('branch_ratio', float),
    ('photon_branch_ratio', float),
    ('photon_branch_ratio_err', float),
    ('beta_branch_ratio', float),
    ('beta_branch_ratio_err', float),
])

gammas_dtype = np.dtype([
    ('from_nuc', int),
    ('to_nuc', int),
    ('parent_nuc', int),
    ('energy', float),
    ('energy_err', float),
    ('photon_intensity', float),
    ('photon_intensity_err', float),
    ('conv_intensity', float),
    ('conv_intensity_err', float),
    ('total_intensity', float),
    ('total_intensity_err', float),
    ('k_conv_e', float),
    ('l_conv_e', float),
    ('m_conv_e', float),
])

alphas_dtype = np.dtype([
    ('from_nuc', int),
    ('to_nuc', int),
    ('energy', float),
    ('intensity', float),
])

betas_dtype = np.dtype([
    ('from_nuc', int),
    ('to_nuc', int),
    ('endpoint_energy', float),
    ('avg_energy', float),
    ('intensity', float),
])

ecbp_dtype = np.dtype([
    ('from_nuc', int),
    ('to_nuc', int),
    ('endpoint_energy', float),
    ('avg_energy', float),
    ('beta_plus_intensity', float),
    ('ec_intensity', float),
    ('k_conv_e', float),
    ('l_conv_e', float),
    ('m_conv_e', float),
])


def parse_level_data(build_dir=""):
    """Builds and returns a list of nuclide decay data.
    Parameters
    ----------
    build_dir : str
        build_nuc_data directory containing ENSDF folder

    Returns
    -------
    level_list_array : np.ndarray
        array of level data
    """
    build_dir = os.path.join(build_dir, 'ENSDF')

    level_list = []
    files = sorted([f for f in glob.glob(os.path.join(build_dir, 'ensdf.*'))])
    for f in files:
        print("    building level data from {0}".format(f))
        level_list = ensdf.levels(f, level_list)

    level_list_array = np.array(level_list, dtype=level_dtype)

    return level_list_array


def parse_decay_data(build_dir=""):
    """Builds and returns a list of nuclide decay data.
    Parameters
    ----------
    build_dir : str
        build_nuc_data directory containing ENSDF folder

    Returns
    -------
    all_decay_array : np.ndarray
        array of decay data
    all_gammas_array : np.ndarray
        array of gamma ray data
    all_alphas_array : np.ndarray
        array of alpha decay data
    all_betas_array : np.ndarray
        array of beta decay data
    all_ecbp_array : np.ndarray
        array of electron capture and beta plus decay data
    """
    build_dir = os.path.join(build_dir, 'ENSDF')

    decay_data = []
    files = sorted([f for f in glob.glob(os.path.join(build_dir, 'ensdf.*'))])
    for f in files:
        print("    parsing decay data from {0}".format(f))
        decay_data = ensdf.decays(f, decay_data)

    all_decays = []
    all_gammas = []
    all_alphas = []
    all_betas = []
    all_ecbp = []
    for item in decay_data:
        all_decays.append(item[:10])
        if len(item[10]) > 0:
            for subitem in item[10]:
                all_gammas.append(tuple(subitem))
        if len(item[11]) > 0:
            for subitem in item[11]:
                all_alphas.append(tuple(subitem))
        if len(item[12]) > 0:
            for subitem in item[12]:
                all_betas.append(tuple(subitem))
        if len(item[13]) > 0:
            for subitem in item[13]:
                all_ecbp.append(tuple(subitem))

    all_decay_array = np.array(all_decays, dtype=decay_dtype)
    all_gammas_array = np.array(all_gammas, dtype=gammas_dtype)
    all_alphas_array = np.array(all_alphas, dtype=alphas_dtype)
    all_betas_array = np.array(all_betas, dtype=betas_dtype)
    all_ecbp_array = np.array(all_ecbp, dtype=ecbp_dtype)

    return all_decay_array, all_gammas_array, all_alphas_array, \
           all_betas_array, all_ecbp_array


def make_decay_half_life_table(nuc_data, build_dir=""):
    """Makes a decay table in the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place ensdf files in.
    """
    # Grab raw level data
    level_list = parse_level_data(build_dir)

    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS)

    # Make a new the table
    if not hasattr(db.root, 'decay'):
        db.createGroup('/', 'decay', 'ENSDF Decay data')

    ll_table = db.createTable('/decay/', 'level_list', level_list,
                              'nuclide [nuc_id], level [keV], half life [s],'
                              'metastable [int]', expectedrows=len(level_list))
    ll_table.flush()

    #now that the level data is in nuc_data we can build the decay data fast
    decay, gammas, alphas, betas, ecbp = parse_decay_data(build_dir)

    decay_table = db.createTable('/decay/', 'decays', decay,
                                 'parent nuclide [nuc_id], daughter nuclide '
                                 '[nuc_id], decay [string], half life [s],'
                                 'half life error [s], branch ratio [frac],'
                                 'photon branch ratio [ratio],'
                                 'photon branch ratio error [ratio],'
                                 'beta branch ratio [ratio],'
                                 'beta branch ratio error [ratio]',
                                 expectedrows=len(decay))
    decay_table.flush()

    gamma_table = db.createTable('/decay/', 'gammas', gammas,
                                 'from_nuc [int], to_nuc [int], primary parent'
                                 'nuc_id [int],'
                                 'Energy [keV], Energy error [keV], '
                                 'photon intensity [ratio], '
                                 'photon intensity error [ratio],'
                                 'conversion e intensity [ratio],'
                                 'conversion e intensity error [ratio],'
                                 'total intensity [ratio],'
                                 'total intensity error [ratio], '
                                 'K conversion electron'
                                 'intensity [ratio], L conversion electron'
                                 'intensity [ratio], M conversion electron'
                                 'intensity [ratio]',
                                 expectedrows=len(gammas))

    gamma_table.flush()

    alphas_table = db.createTable('/decay/', 'alphas', alphas,
                                  'from_nuc [int], to_nuc [int]'
                                  'Energy [keV], Intensity [ratio],',
                                  expectedrows=len(alphas))
    alphas_table.flush()

    betas_table = db.createTable('/decay/', 'betas', betas,
                                 'from_nuc [int], to_nuc [int],'
                                 'Endpoint Energy [keV], Average Energy [keV],'
                                 'Intensity [ratio]'
                                 , expectedrows=len(betas))

    betas_table.flush()

    ecbp_table = db.createTable('/decay/', 'ecbp', ecbp,
                                'from_nuc [int], to_nuc [int],'
                                'Endpoint Energy [keV], Average Energy [keV],'
                                'B+ Intensity [ratio], '
                                'Electron Capture Intensity [ratio],'
                                'K conversion'
                                'electron intensity [ratio], L conversion'
                                'electron intensity [ratio], M conversion'
                                'electron intensity [ratio]',
                                expectedrows=len(ecbp))
    ecbp_table.flush()

    # Close the hdf5 file
    db.close()


def make_decay(args):
    """Controller function for adding decay data."""
    nuc_data, build_dir = args.nuc_data, args.build_dir

    with tb.openFile(nuc_data, 'r') as f:
        if hasattr(f.root, 'decay'):
            print("skipping ENSDF decay data table creation; already exists.")
            return

            # grab the decay data
    print("Grabbing the ENSDF decay data from IAEA")
    grab_ensdf_decay(build_dir)

    # Make atomic mass table once we have the array
    print("Making decay data table.")
    make_decay_half_life_table(nuc_data, build_dir)

