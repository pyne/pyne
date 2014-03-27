"""This module provides a way to grab and store raw data for radioactive decay."""
from __future__ import print_function
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
    ensdf_zip = ['ensdf_131009_099.zip', 'ensdf_131009_199.zip', 'ensdf_131009_294.zip',]

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


half_life_dtype = np.dtype([
    ('from_nuc', int),
    ('to_nuc', int),
    ('level', float),
    ('half_life', float),
    ('decay_const', float),
    ('branch_ratio', float),
    ])

level_dtype = np.dtype([
    ('nuc_id', int),
    ('level', float),
    ('half_life', float),
    ('metastable', int),
    ])

decay_dtype = np.dtype([
    ('parent', int),
    ('child', int),
    ('decay', 'S5'),
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


def parse_decay(build_dir=""):
    """Builds and returns a list of nuclide decay data.
    Parameters
    ----------
    build_dir : str
        build_nuc_data directory containing ENSDF folder

    Returns
    -------
    half_life_data_array : np.ndarray
        array of half life data
    level_list_array : np.ndarray
        array of level data
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

    half_life_data = []
    decay_data = []
    level_list = []
    lmap = dict()
    lcount = 0
    files = sorted([f for f in glob.glob(os.path.join(build_dir, 'ensdf.*'))])
    for f in files:
        print("    parsing decay data from {0}".format(f))
        half_life_data += ensdf.half_life(f)
        level_list, decay_data, lmap, lcount = \
            ensdf.decays(f, level_list, decay_data, lmap, lcount)

    ln2 = np.log(2.0)
    half_life_data = [(fn, tn, lvl, hl, ln2 / hl, br)
                      for fn, lvl, tn, hl, br in half_life_data]
    half_life_data = set(half_life_data)
    half_life_data = sorted(half_life_data, key=lambda x: (x[0], x[1]))

    half_life_data_array = np.array(half_life_data, dtype=half_life_dtype)
    #da, mask = np.unique(decay_array, return_index=True)
    #mask.sort()
    #decay_array = decay_array[mask]
    level_list_array = np.array(level_list, dtype=level_dtype)
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

    return half_life_data_array, level_list_array, all_decay_array, \
           all_gammas_array, all_alphas_array, all_betas_array, all_ecbp_array


def make_decay_half_life_table(nuc_data, build_dir=""):
    """Makes a decay table in the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place ensdf files in.
    """
    # Grab raw data
    half_life, level_list, decay, \
    gammas, alphas, betas, ecbp = parse_decay(build_dir)

    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS)

    # Make a new the table
    if not hasattr(db.root, 'decay'):
        db.createGroup('/', 'decay', 'ENSDF Decay data')

    decaytable = db.createTable('/decay/', 'half_life',
                                np.empty(0, dtype=half_life_dtype),
                                'Decay Energy level [MeV], half_life [s],'
                                'decay_const [1/s], branch_ratio [frac]',
                                expectedrows=len(half_life))
    decaytable.append(half_life)

    # Ensure that data was written to table
    decaytable.flush()

    ll_table = db.createTable('/decay/', 'level_list', level_list,
                              'nuclide [nuc_id], level [keV], half life [s],'
                              'metastable [int]', expectedrows=len(level_list))
    ll_table.flush()

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
                                 'Energy [keV], Energy error [keV], '
                                 'photon intensity [ratio], '
                                 'photon intensity error [ratio],'
                                 'conversion e intensity [ratio],'
                                 'conversion e intensity error [ratio],'
                                 'total intensity [ratio],'
                                 'total intensity error [ratio], '
                                 'from_nuc [int], to_nuc [int], primary parent'
                                 'nuc_id [int], K conversion electron'
                                 'intensity [ratio], L conversion electron'
                                 'intensity [ratio], M conversion electron'
                                 'intensity [ratio]',
                                 expectedrows=len(gammas))

    gamma_table.flush()

    alphas_table = db.createTable('/decay/', 'alphas', alphas,
                                  'Energy [keV], Intensity [ratio],'
                                  'from_nuc [int], to_nuc [int]',
                                  expectedrows=len(alphas))
    alphas_table.flush()

    betas_table = db.createTable('/decay/', 'betas', betas,
                                 'Endpoint Energy [keV], Average Energy [keV],'
                                 'Intensity [ratio], from_nuc [int],'
                                 'to_nuc [int]', expectedrows=len(betas))

    betas_table.flush()

    ecbp_table = db.createTable('/decay/', 'ecbp', ecbp,
                                'Endpoint Energy [keV], Average Energy [keV],'
                                'B+ Intensity [ratio], '
                                'Electron Capture Intensity [ratio],'
                                'from_nuc [int], to_nuc [int], K conversion'
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

