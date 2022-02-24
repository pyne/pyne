"""This module provides a way to grab and store raw data for neutron scattering 
lengths.  This data comes from Neutron News, Vol. 3, No. 3, 1992, pp. 29-37 via 
a NIST webpage (http://www.ncnr.nist.gov/resources/n-lengths/list.html).  Please
contact Alan Munter, <alan.munter@nist.gov> for more information."""
from __future__ import print_function
import os
import re
import shutil
from pyne.utils import QA_warn

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2

import numpy as np
import tables as tb

from .. import nucname
from .api import BASIC_FILTERS

QA_warn(__name__)


def grab_scattering_lengths(build_dir="", file_out="scattering_lengths.html"):
    """Grabs the scattering cross-section lengths for neutrons from the NIST website
    or locally from this module."""
    build_filename = os.path.join(build_dir, file_out)
    local_filename = os.path.join(os.path.split(__file__)[0], file_out)

    if os.path.exists(local_filename):
        shutil.copy(local_filename, build_filename)
        return

    nist = urllib2.urlopen("http://www.ncnr.nist.gov/resources/n-lengths/list.html")
    with open(build_filename, "w") as f:
        f.write(nist.read())


def nist_num(nist_data):
    """Converts a NIST style data point to a point.

    Parameters
    ----------
    nist_data : str
        A nist data point.

    Returns
    -------
    d : float or complex
        a data point.
    """
    nd = nist_data
    while ("(" in nd) or (")" in nd):
        nd_pre = nd.partition("(")
        nd_post = nd.partition(")")
        nd = nd_pre[0] + nd_post[2]

    if (nd == "---") or (nd == ""):
        nd = "0.0"

    nd = nd.replace("<i>i</i>", "j")
    d = eval(nd)
    return d


sl_dtype = np.dtype(
    [
        ("nuc", int),
        ("b_coherent", np.complex128),
        ("b_incoherent", np.complex128),
        ("xs_coherent", float),
        ("xs_incoherent", float),
        ("xs", float),
    ]
)

scat_len_data = "[ aeE()<>i/.+\d-]+?"
scat_len_space = "[ \t]+"
scat_len_pattern = "<td>{space}(?P<iso>[A-Za-z\d]+){space}<td>{space}(?P<conc>{data}){space}<td>{space}(?P<b_coherent>{data}){space}<td>{space}(?P<b_incoherent>{data}){space}<td>{space}(?P<xs_coherent>{data}){space}<td>{space}(?P<xs_incoherent>{data}){space}<td>{space}(?P<xs>{data}){space}<td>{space}(?P<xs_a>{data}){space}<tr>".format(
    data=scat_len_data, space=scat_len_space
)


def parse_scattering_lengths(build_dir):
    """Converts to scattering lenth data to a numpy array."""
    build_filename = os.path.join(build_dir, "scattering_lengths.html")

    # Read in cinder data file
    with open(build_filename, "r") as f:
        raw_data = f.read()

    sl_data = []

    # Iterate over all isotopes in the table
    for m in re.finditer(scat_len_pattern, raw_data):
        md = m.groupdict()

        slrow = (
            nucname.id(md["iso"]),
            nist_num(md["b_coherent"]) * (1e-13),
            nist_num(md["b_incoherent"]) * (1e-13),
            nist_num(md["xs_coherent"]),
            nist_num(md["xs_incoherent"]),
            nist_num(md["xs"]),
        )

        sl_data.append(slrow)

    sl_array = np.array(sl_data, dtype=sl_dtype)
    return sl_array


def make_scattering_lengths_table(nuc_data, build_dir=""):
    """Adds the neutron sacttering lengths to the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place html files in.
    """
    sl_array = parse_scattering_lengths(build_dir)

    # Open the HDF5 File
    db = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Ensure that the appropriate file structure is present
    if not hasattr(db.root, "neutron"):
        # Create neutron group
        neutron_group = db.create_group("/", "neutron", "Neutron Data")

    # Init the neutron fission product info table
    sl_table = db.create_table(
        "/neutron/",
        "scattering_lengths",
        np.empty(0, dtype=sl_dtype),
        "Neutron Scattering Lengths, b [cm], sigma [barns]",
        expectedrows=len(sl_array),
    )
    sl_table.append(sl_array)

    # Write the table
    sl_table.flush()

    # Close the hdf5 file
    db.close()


def make_scattering_lengths(args):
    """Controller function for adding scattering lengths."""
    nuc_data, build_dir = args.nuc_data, args.build_dir

    # Check that the table exists
    with tb.open_file(nuc_data, "a", filters=BASIC_FILTERS) as f:
        if hasattr(f.root, "neutron") and hasattr(f.root.neutron, "scattering_lengths"):
            print("skipping scattering lengths data table creation; already exists.")
            return

    # Grab the raw data
    print("Grabbing the scattering length data.")
    grab_scattering_lengths(build_dir)

    # Make scatering table once we have the data
    print("Making neutron scattering length table.")
    make_scattering_lengths_table(nuc_data, build_dir)
