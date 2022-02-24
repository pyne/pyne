"""Module handles parsing EAF formatted cross section files and adding
the data to PyNE's HDF5 storage.  The data here is autonatically grabbed from
the IAEA. 
"""
from __future__ import print_function
import re
import os
from pyne.utils import QA_warn

try:
    import urllib.request as urllib
except ImportError:
    import urllib
from gzip import GzipFile

import numpy as np
import tables as tb

from .. import nucname
from .api import BASIC_FILTERS

QA_warn(__name__)


def grab_eaf_data(build_dir=""):
    """Grabs the EAF activation data files
    if not already present.

    Parameters
    ----------
    build_dir : str
        Major directory to place EAF data file(s) in. 'EAF/' will be appended.
    """
    # Add EAF to build_dir
    build_dir = os.path.join(build_dir, "EAF")
    try:
        os.makedirs(build_dir)
        print("{0} created".format(build_dir))
    except OSError:
        pass

    # Grab ENSDF files and unzip them.
    # This link was taken from 'http://www-nds.iaea.org/fendl/fen-activation.htm'
    iaea_url = "http://www-nds.iaea.org/fendl2/activation/processed/vitj_e/libout/fendlg-2.0_175-gz"
    cf_base_url = "https://github.com/pyne/data/raw/master/"
    eaf_gzip = "fendlg-2.0_175-gz"

    fpath = os.path.join(build_dir, eaf_gzip)
    if eaf_gzip not in os.listdir(build_dir):
        print("  grabbing {0} and placing it in {1}".format(eaf_gzip, fpath))
        try:
            urllib.urlretrieve(iaea_url, fpath)
        except (OSError, IOError):
            open(fpath, "a").close()  # touch the file

        if os.path.getsize(fpath) < 3215713:
            print("  could not get {0} from IAEA; trying S3 mirror".format(eaf_gzip))
            os.remove(fpath)
            try:
                urllib.urlretrieve(cf_base_url + eaf_gzip, fpath)
            except (OSError, IOError):
                open(fpath, "a").close()  # touch the file
            if os.path.getsize(fpath) < 3215713:
                print("  could not get {0} from S3 mirror".format(eaf_gzip))
                return False

    # Write contents of single-file gzip archive to a new file
    try:
        gf = GzipFile(fpath)
        ofile = os.path.join(build_dir, "fendlg-2.0_175")
        with open(ofile, "w") as fw:
            for line in gf:
                fw.write(line.decode("us-ascii"))
    finally:
        gf.close()

    return True


# numpy array row storage information for EAF data
eaf_dtype = np.dtype(
    [
        ("nuc_zz", int),
        ("rxnum", "S7"),
        ("rxstr", "S7"),
        ("daughter", "S7"),
        ("xs", float, (175,)),
    ]
)

# Regular expression for parsing an individual set of EAF data.
# Includes some groupnames that are currently unused.
eaf_info_pattern = (
    "(?P<iso>\d{5,7})\s*"
    + "(?P<rxnum>\d{2,4})\s*"
    + "(?P<ngrps>\d{1,3})\s*"
    + "(?P<parent>[a-zA-Z]{1,2}\s{0,3}\d{1,3}[M ][12 ])"
    + "(?P<rxstr>\(N,[\w\s]{3}\))"
    + "(?P<daugh>[a-zA-Z.]{1,2}\s{0,3}\d{1,3}[MG]{0,1}\d{0,1})"
    + "(.*?)"
)
eaf_bin_pattern = "(?P<xs>(\d\.\d{5}E[-+]\d{2}\s*){1,175})"


def parse_eaf_xs(build_file):
    """Create numpy array by parsing EAF data
    using regular expressions

    Parameters
    ----------
    build_file : str
        Path where EAF data is stored.

    Returns
    -------
    eaf_array : numpy array
        Numpy array with a row for each isotope+reaction combination
        found in the EAF data.

    """

    with open(build_file, "r") as f:
        raw_data = f.read()

    eaf_data = list()
    eaf_pattern = eaf_info_pattern + eaf_bin_pattern

    # Iterate over all iso/rx combinations in file
    for m in re.finditer(eaf_pattern, raw_data, re.DOTALL):
        md = m.groupdict()

        xs_list = [float(x) for x in md["xs"].split()]
        xs_list += (175 - len(xs_list)) * [0.0]

        # Store information in new row of array.
        eafrow = (nucname.id(md["iso"]), md["rxnum"], md["rxstr"], md["daugh"], xs_list)

        eaf_data.append(eafrow)

    eaf_array = np.array(eaf_data, dtype=eaf_dtype)

    print("Read in {0} sets of EAF data.".format(len(eaf_array)))

    return eaf_array


def make_eaf_table(nuc_data, build_path=""):
    """Function for adding EAF group and table to HDF5 storage.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_path : str
        Directory where EAF data is located.

    """

    print("Grabbing the EAF activation data.")
    eaf_array = parse_eaf_xs(build_path)

    # Open the HDF5 file
    db = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Ensure that the appropriate file structure is present
    if not hasattr(db.root, "neutron"):
        neutron_group = db.create_group("/", "neutron", "Neutron Interaction Data")

    # Create eaf_xs group
    if not hasattr(db.root.neutron, "eaf_xs"):
        eaf_group = db.create_group(
            "/neutron", "eaf_xs", "EAF 175-Group Neutron Activation Cross Section Data"
        )

    eaf_table = db.create_table(
        "/neutron/eaf_xs",
        "eaf_xs",
        np.empty(0, dtype=eaf_dtype),
        "EAF Activation Cross Section Data [barns]",
        expectedrows=len(eaf_array),
    )

    # Put eaf_array in the table
    eaf_table.append(eaf_array)

    # Write the table
    eaf_table.flush()

    # Add group structure by calling placeholder function to get boundaries
    db.create_array(
        "/neutron/eaf_xs", "E_g", _get_eaf_groups(), "Neutron energy group bounds [MeV]"
    )

    # Close the HDF5 file
    db.close()


def _get_eaf_groups():
    """Function for retrieving EAF group structure.

    Note, energies are entered high to low in the HDF5 storage.

    This is used as a placeholder for a future PyNE implementation of a database
    storing group structures.  This hypothetical database would return the
    group structure like this method does.

    Returns
    ----------
    eaf_E_g_array : 1D list
        List of energy group boundaries from high to low.

    """

    eaf_E_g_array = [
        0.0,
        1.0000e-07,
        4.1399e-07,
        5.3158e-07,
        6.8256e-07,
        8.7643e-07,
        1.1254e-06,
        1.4450e-06,
        1.8554e-06,
        2.3824e-06,
        3.0590e-06,
        3.9279e-06,
        5.0435e-06,
        6.4760e-06,
        8.3153e-06,
        1.0677e-05,
        1.3710e-05,
        1.7604e-05,
        2.2603e-05,
        2.9023e-05,
        3.7267e-05,
        4.7851e-05,
        6.1442e-05,
        7.8893e-05,
        1.0130e-04,
        1.3007e-04,
        1.6702e-04,
        2.1445e-04,
        2.7536e-04,
        3.5358e-04,
        4.5400e-04,
        5.8295e-04,
        7.4852e-04,
        9.6112e-04,
        1.2341e-03,
        1.5846e-03,
        2.0347e-03,
        2.2487e-03,
        2.4852e-03,
        2.6126e-03,
        2.7465e-03,
        3.0354e-03,
        3.3546e-03,
        3.7074e-03,
        4.3074e-03,
        5.5308e-03,
        7.1017e-03,
        9.1188e-03,
        1.0595e-02,
        1.1709e-02,
        1.5034e-02,
        1.9305e-02,
        2.1875e-02,
        2.3579e-02,
        2.4176e-02,
        2.4788e-02,
        2.6058e-02,
        2.7000e-02,
        2.8501e-02,
        3.1828e-02,
        3.4307e-02,
        4.0868e-02,
        4.6309e-02,
        5.2475e-02,
        5.6562e-02,
        6.7380e-02,
        7.2025e-02,
        7.9499e-02,
        8.2503e-02,
        8.6517e-02,
        9.8037e-02,
        1.1109e-01,
        1.1679e-01,
        1.2277e-01,
        1.2907e-01,
        1.3569e-01,
        1.4264e-01,
        1.4996e-01,
        1.5764e-01,
        1.6573e-01,
        1.7422e-01,
        1.8316e-01,
        1.9255e-01,
        2.0242e-01,
        2.1280e-01,
        2.2371e-01,
        2.3518e-01,
        2.4724e-01,
        2.7324e-01,
        2.8725e-01,
        2.9452e-01,
        2.9721e-01,
        2.9849e-01,
        3.0197e-01,
        3.3373e-01,
        3.6883e-01,
        3.8774e-01,
        4.0762e-01,
        4.5049e-01,
        4.9787e-01,
        5.2340e-01,
        5.5023e-01,
        5.7844e-01,
        6.0810e-01,
        6.3928e-01,
        6.7206e-01,
        7.0651e-01,
        7.4274e-01,
        7.8082e-01,
        8.2085e-01,
        8.6294e-01,
        9.0718e-01,
        9.6167e-01,
        1.0026e00,
        1.1080e00,
        1.1648e00,
        1.2246e00,
        1.2874e00,
        1.3534e00,
        1.4227e00,
        1.4957e00,
        1.5724e00,
        1.6530e00,
        1.7377e00,
        1.8268e00,
        1.9205e00,
        2.0190e00,
        2.1225e00,
        2.2313e00,
        2.3069e00,
        2.3457e00,
        2.3653e00,
        2.3851e00,
        2.4660e00,
        2.5924e00,
        2.7253e00,
        2.8651e00,
        3.0119e00,
        3.1664e00,
        3.3287e00,
        3.6788e00,
        4.0657e00,
        4.4933e00,
        4.7237e00,
        4.9659e00,
        5.2205e00,
        5.4881e00,
        5.7695e00,
        6.0653e00,
        6.3763e00,
        6.5924e00,
        6.7032e00,
        7.0469e00,
        7.4082e00,
        7.7880e00,
        8.1873e00,
        8.6071e00,
        9.0484e00,
        9.5123e00,
        1.0000e01,
        1.0513e01,
        1.1052e01,
        1.1618e01,
        1.2214e01,
        1.2523e01,
        1.2840e01,
        1.3499e01,
        1.3840e01,
        1.4191e01,
        1.4550e01,
        1.4918e01,
        1.5683e01,
        1.6487e01,
        1.6905e01,
        1.7333e01,
        1.9640e01,
    ]

    eaf_E_g_array.reverse()

    return eaf_E_g_array


def make_eaf(args):
    """Controller function for adding cross section data from EAF format file."""
    nuc_data, build_dir, datapath = args.nuc_data, args.build_dir, args.datapath

    # Check if the table already exists
    with tb.open_file(nuc_data, "a", filters=BASIC_FILTERS) as f:
        if hasattr(f.root, "neutron") and hasattr(f.root.neutron, "eaf_xs"):
            print("skipping EAF activation data table creation; already exists.")
            return

    # grab the EAF data
    print("Grabbing the EAF activation data from IAEA")
    grabbed = grab_eaf_data(build_dir)

    if not grabbed:
        return

    build_filename = os.path.join(build_dir, "EAF/fendlg-2.0_175")

    if os.path.exists(build_filename):
        build_path = build_filename
    else:
        return

    #
    print("Making EAF activation data table.")
    make_eaf_table(nuc_data, build_path)
