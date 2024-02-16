"""This module provides a way to grab and store raw data for radioactive decay."""
from __future__ import print_function, division
import os
import glob
import shutil
from zipfile import ZipFile
from pyne.utils import QA_warn

try:
    import urllib.request as urllib
except ImportError:
    import urllib2 as urllib


import numpy as np
import tables as tb

from pyne import ensdf
from pyne.dbgen.api import BASIC_FILTERS

QA_warn(__name__)


def _readpoint(line, dstart, dlen):
    data = ensdf._getvalue(line[dstart : dstart + dlen])
    error = ensdf._getvalue(line[dstart + dlen : dstart + dlen + 2])
    return data, error


def _read_variablepoint(line, dstart, dlen):
    sub = line[dstart : dstart + dlen + 2].split()
    data = None
    error = None
    if len(sub) == 2:
        data = ensdf._getvalue(sub[0])
        error = ensdf._getvalue(sub[1])
    return data, error


def grab_atomic_data(build_dir=""):
    medfile = os.path.join(build_dir, "mednew.dat")
    local = os.path.join(os.path.dirname(__file__), "mednew.dat")
    if os.path.isfile(medfile):
        return
    # try to download first
    # nndc url seems to be down
    # url = 'http://www.nndc.bnl.gov/nndcscr/ensdf_pgm/analysis/radlst/mednew.dat'
    url = "https://www-nds.iaea.org/workshops/smr1939/Codes/ENSDF_Codes/mswindows/radlst/mednew.dat"
    try:
        conn = urllib.urlopen(url)
        with open(medfile, "wb") as f:
            f.write(conn.read())
        return
    except Exception:
        pass
    # use local copy if we can't download
    shutil.copy(local, medfile)


def parse_atomic_data(build_dir=""):
    i = 0
    j = 0
    medfile = os.path.join(build_dir, "mednew.dat")
    dat = np.zeros((103,), atomic_dtype)
    with open(medfile, "r") as f:
        lines = f.readlines()
    for line in lines:
        if (-1) ** i == 1:
            Z = int(line[0:3])
            k_shell_fluor, k_shell_fluor_error = _readpoint(line, 9, 6)
            l_shell_fluor, l_shell_fluor_error = _readpoint(line, 18, 6)
            # Probability of creating L-shell vacancy by filling K-shell
            prob, prob_error = _readpoint(line, 27, 6)
            k_shell_be, k_shell_be_err = _readpoint(line, 36, 8)
            li_shell_be, li_shell_be_err = _readpoint(line, 47, 8)
            mi_shell_be, mi_shell_be_err = _readpoint(line, 58, 8)
            ni_shell_be, ni_shell_be_err = _readpoint(line, 69, 8)
        else:
            Kb_to_Ka, Kb_to_Ka_err = _read_variablepoint(line, 9, 7)
            Ka2_to_Ka1, Ka2_to_Ka1_err = _read_variablepoint(line, 19, 7)
            L_auger = ensdf._getvalue(line[29:36])
            K_auger = ensdf._getvalue(line[36:42])
            Ka1_X_ray_en, Ka1_X_ray_en_err = _readpoint(line, 43, 8)
            Ka2_X_ray_en, Ka2_X_ray_en_err = _readpoint(line, 54, 7)
            Kb_X_ray_en = ensdf._getvalue(line[65:69])
            L_X_ray_en = ensdf._getvalue(line[70:76])
            dat[j] = (
                Z,
                k_shell_fluor,
                k_shell_fluor_error,
                l_shell_fluor,
                l_shell_fluor_error,
                prob,
                k_shell_be,
                k_shell_be_err,
                li_shell_be,
                li_shell_be_err,
                mi_shell_be,
                mi_shell_be_err,
                ni_shell_be,
                ni_shell_be_err,
                Kb_to_Ka,
                Kb_to_Ka_err,
                Ka2_to_Ka1,
                Ka2_to_Ka1_err,
                L_auger,
                K_auger,
                Ka1_X_ray_en,
                Ka1_X_ray_en_err,
                Ka2_X_ray_en,
                Ka2_X_ray_en_err,
                Kb_X_ray_en,
                L_X_ray_en,
            )
            j += 1
        i += 1
    return dat


def grab_ensdf_decay(build_dir=""):
    """
    Grabs the ENSDF decay data files
    if not already present.

    Parameters
    ----------
    build_dir : str
        Major directory to place html files in. 'ENSDF/' will be appended.
    """
    # Add ENSDF to build_dir
    build_dir = os.path.join(build_dir, "ENSDF")
    try:
        os.makedirs(build_dir)
    except OSError:
        pass

    # Grab ENSDF files and unzip them.
    iaea_base_url = "http://www.nndc.bnl.gov/ensarchivals/distributions/dist19/"

    cf_base_url = "https://github.com/pyne/data/raw/master/"
    ensdf_zip = [
        "ensdf_191004_099.zip",
        "ensdf_191004_199.zip",
        "ensdf_191004_300.zip",
    ]

    for f in ensdf_zip:
        fpath = os.path.join(build_dir, f)
        if f not in os.listdir(build_dir):
            print("  grabbing {0} and placing it in {1}".format(f, fpath))
            conn = urllib.urlopen(iaea_base_url + f)
            with open(fpath, "wb") as f:
                f.write(conn.read())

            if os.path.getsize(fpath) < 1048576:
                print("  could not get {0} from NNDC; trying mirror".format(f))
                os.remove(fpath)
                conn = urllib.urlopen(cf_base_url + f)
                with open(fpath, "wb") as f:
                    f.write(conn.read())

        # not using ZipFile context manager (with statement for Python 2.6)
        try:
            zf = ZipFile(fpath)
            for name in zf.namelist():
                if not os.path.exists(os.path.join(build_dir, name)):
                    print("    extracting {0} from {1}".format(name, fpath))
                    zf.extract(name, build_dir)
        finally:
            zf.close()


level_dtype = np.dtype(
    [
        ("nuc_id", int),
        ("rx_id", np.uint32),
        ("half_life", float),
        ("level", float),
        ("branch_ratio", float),
        ("metastable", int),
        ("special", "S1"),
    ]
)

decay_dtype = np.dtype(
    [
        ("parent", int),
        ("child", int),
        ("decay", np.uint32),
        ("half_life", float),
        ("half_life_error", float),
        ("branch_ratio", float),
        ("branch_ratio_error", float),
        ("photon_branch_ratio", float),
        ("photon_branch_ratio_err", float),
        ("beta_branch_ratio", float),
        ("beta_branch_ratio_err", float),
    ]
)

gammas_dtype = np.dtype(
    [
        ("from_nuc", int),
        ("to_nuc", int),
        ("parent_nuc", int),
        ("child_nuc", int),
        ("energy", float),
        ("energy_err", float),
        ("photon_intensity", float),
        ("photon_intensity_err", float),
        ("conv_intensity", float),
        ("conv_intensity_err", float),
        ("total_intensity", float),
        ("total_intensity_err", float),
        ("k_conv_e", float),
        ("l_conv_e", float),
        ("m_conv_e", float),
    ]
)

alphas_dtype = np.dtype(
    [
        ("from_nuc", int),
        ("to_nuc", int),
        ("energy", float),
        ("intensity", float),
    ]
)

betas_dtype = np.dtype(
    [
        ("from_nuc", int),
        ("to_nuc", int),
        ("endpoint_energy", float),
        ("avg_energy", float),
        ("intensity", float),
    ]
)

ecbp_dtype = np.dtype(
    [
        ("from_nuc", int),
        ("to_nuc", int),
        ("endpoint_energy", float),
        ("avg_energy", float),
        ("beta_plus_intensity", float),
        ("ec_intensity", float),
        ("k_conv_e", float),
        ("l_conv_e", float),
        ("m_conv_e", float),
    ]
)


atomic_dtype = np.dtype(
    [
        ("z", int),
        ("k_shell_fluor", float),
        ("k_shell_fluor_error", float),
        ("l_shell_fluor", float),
        ("l_shell_fluor_error", float),
        ("prob", float),
        ("k_shell_be", float),
        ("k_shell_be_err", float),
        ("li_shell_be", float),
        ("li_shell_be_err", float),
        ("mi_shell_be", float),
        ("mi_shell_be_err", float),
        ("ni_shell_be", float),
        ("ni_shell_be_err", float),
        ("kb_to_ka", float),
        ("kb_to_ka_err", float),
        ("ka2_to_ka1", float),
        ("ka2_to_ka1_err", float),
        ("l_auger", float),
        ("k_auger", float),
        ("ka1_x_ray_en", float),
        ("ka1_x_ray_en_err", float),
        ("ka2_x_ray_en", float),
        ("ka2_x_ray_en_err", float),
        ("kb_x_ray_en", float),
        ("l_x_ray_en", float),
    ]
)


def parse_level_data(build_dir=""):
    """
    Builds and returns a list of nuclide decay data.

    Parameters
    ----------
    build_dir : str
        build_nuc_data directory containing ENSDF folder

    Returns
    -------
    level_list_array : np.ndarray
        array of level data
    """
    build_dir = os.path.join(build_dir, "ENSDF")

    level_list = []
    files = sorted([f for f in glob.glob(os.path.join(build_dir, "ensdf.*"))])
    for f in files:
        print("    building level data from {0}".format(f))
        level_list = ensdf.levels(f, level_list)

    level_list_array = np.array(level_list, dtype=level_dtype)

    return level_list_array


def parse_decay_data(build_dir=""):
    """
    Builds and returns a list of nuclide decay data.

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
    build_dir = os.path.join(build_dir, "ENSDF")

    decay_data = []
    files = sorted([f for f in glob.glob(os.path.join(build_dir, "ensdf.*"))])
    for f in files:
        print("    parsing decay data from {0}".format(f))
        decay_data = ensdf.decays(f, decay_data)

    all_decays = []
    all_gammas = []
    all_alphas = []
    all_betas = []
    all_ecbp = []
    for item in decay_data:
        all_decays.append(item[:11])
        if len(item[11]) > 0:
            for subitem in item[11]:
                all_gammas.append(tuple(subitem))
        if len(item[12]) > 0:
            for subitem in item[12]:
                all_alphas.append(tuple(subitem))
        if len(item[13]) > 0:
            for subitem in item[13]:
                all_betas.append(tuple(subitem))
        if len(item[14]) > 0:
            for subitem in item[14]:
                all_ecbp.append(tuple(subitem))

    all_decay_array = np.array(all_decays, dtype=decay_dtype)
    all_gammas_array = np.array(all_gammas, dtype=gammas_dtype)
    all_alphas_array = np.array(all_alphas, dtype=alphas_dtype)
    all_betas_array = np.array(all_betas, dtype=betas_dtype)
    all_ecbp_array = np.array(all_ecbp, dtype=ecbp_dtype)

    return (
        all_decay_array,
        all_gammas_array,
        all_alphas_array,
        all_betas_array,
        all_ecbp_array,
    )


def make_atomic_decay_table(nuc_data, build_dir=""):
    """Makes atomic decay table in the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place xray data file in.
    """
    xrd = parse_atomic_data(build_dir)

    db = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Make a new the table
    if not hasattr(db.root, "decay"):
        db.create_group("/", "decay", "ENSDF Decay data")

    atomic_table = db.create_table(
        "/decay/",
        "atomic",
        xrd,
        "z"
        "k_shell_fluor"
        "k_shell_fluor_error"
        "l_shell_fluor"
        "l_shell_fluor_error"
        "prob"
        "k_shell_be"
        "k_shell_be_err"
        "li_shell_be"
        "li_shell_be_err"
        "mi_shell_be"
        "mi_shell_be_err"
        "ni_shell_be"
        "ni_shell_be_err"
        "kb_to_ka"
        "kb_to_ka_err"
        "ka2_to_ka1"
        "ka2_to_ka1_err"
        "k_auger"
        "k_auger"
        "ka1_x_ray_en"
        "ka1_x_ray_en_err"
        "ka2_x_ray_en"
        "ka2_x_ray_en_err"
        "kb_x_ray_en"
        "l_x_ray_en",
        expectedrows=103,
    )
    atomic_table.flush()
    db.close()


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
    db = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Make a new the table
    if not hasattr(db.root, "decay"):
        db.create_group("/", "decay", "ENSDF Decay data")

    ll_table = db.create_table(
        "/decay/",
        "level_list",
        level_list,
        "nuclide [nuc_id], level [keV], half life [s]," "metastable [int]",
        expectedrows=len(level_list),
    )
    ll_table.flush()

    # now that the level data is in nuc_data we can build the decay data fast
    decay, gammas, alphas, betas, ecbp = parse_decay_data(build_dir)

    decay_table = db.create_table(
        "/decay/",
        "decays",
        decay,
        "parent nuclide [nuc_id], daughter nuclide "
        "[nuc_id], decay [string], half life [s],"
        "half life error [s], branch ratio [frac],"
        "branch ratio error [frac],"
        "photon branch ratio [ratio],"
        "photon branch ratio error [ratio],"
        "beta branch ratio [ratio],"
        "beta branch ratio error [ratio]",
        expectedrows=len(decay),
    )
    decay_table.flush()

    gamma_table = db.create_table(
        "/decay/",
        "gammas",
        gammas,
        "from_nuc [int], to_nuc [int], primary parent"
        "nuc_id [int], child nuc_id [int]"
        "Energy [keV], Energy error [keV], "
        "photon intensity [ratio], "
        "photon intensity error [ratio],"
        "conversion e intensity [ratio],"
        "conversion e intensity error [ratio],"
        "total intensity [ratio],"
        "total intensity error [ratio], "
        "K conversion electron"
        "intensity [ratio], L conversion electron"
        "intensity [ratio], M conversion electron"
        "intensity [ratio]",
        expectedrows=len(gammas),
    )

    gamma_table.flush()

    alphas_table = db.create_table(
        "/decay/",
        "alphas",
        alphas,
        "from_nuc [int], to_nuc [int]" "Energy [keV], Intensity [ratio],",
        expectedrows=len(alphas),
    )
    alphas_table.flush()

    betas_table = db.create_table(
        "/decay/",
        "betas",
        betas,
        "from_nuc [int], to_nuc [int],"
        "Endpoint Energy [keV], Average Energy [keV],"
        "Intensity [ratio]",
        expectedrows=len(betas),
    )

    betas_table.flush()

    ecbp_table = db.create_table(
        "/decay/",
        "ecbp",
        ecbp,
        "from_nuc [int], to_nuc [int],"
        "Endpoint Energy [keV], Average Energy [keV],"
        "B+ Intensity [ratio], "
        "Electron Capture Intensity [ratio],"
        "K conversion"
        "electron intensity [ratio], L conversion"
        "electron intensity [ratio], M conversion"
        "electron intensity [ratio]",
        expectedrows=len(ecbp),
    )
    ecbp_table.flush()

    # Close the hdf5 file
    db.close()


def make_decay(args):
    """Controller function for adding decay data."""
    nuc_data, build_dir = args.nuc_data, args.build_dir

    with tb.open_file(nuc_data, "r") as f:
        if hasattr(f.root, "decay") and hasattr(f.root.decay, "ecbp"):
            print("skipping ENSDF decay data table creation; already exists.")
            return

            # grab the decay data
    print("Grabbing the ENSDF decay data from NNDC")
    grab_ensdf_decay(build_dir)

    # Make atomic mass table once we have the array
    print("Making decay data table.")
    make_decay_half_life_table(nuc_data, build_dir)

    print("Grabbing Atomic data")
    grab_atomic_data(build_dir)

    print("Making atomic decay data table")
    make_atomic_decay_table(nuc_data, build_dir)
