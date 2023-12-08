"""This module provides a way to grab and store raw data for fission product yeilds
from the WIMSD library at the IAEA. For more information, please visit their website:
https://www-nds.iaea.org/wimsd/index.html or 
https://www-nds.iaea.org/wimsd/fpyield.htm. Please contact the NDS at  
nds.contact-point@iaea.org with questions about the data itself.

The copyright for the data parsed here is held by the IAEA and is made available 
under the following conditions:

**Disclaimer:** Distributed data products contain consensus values of physical 
constants. However, neither the network centre nor the IAEA guarantees the 
accuracy of such data products or their suitability for particular applied 
scientific purposes.

**Copyright:**  One may use or reproduce data and information from this site with 
an appropriate acknowledgement to the source of data. One may not charge any 
subsequent fee for these data.

"""
from __future__ import print_function
import os
import re
import sys
import shutil
from pyne.utils import QA_warn

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2
if sys.version_info[0] >= 3:
    from html.parser import HTMLParser
else:
    from HTMLParser import HTMLParser

import numpy as np
import tables as tb

from pyne import nucname
from pyne.dbgen.api import BASIC_FILTERS

QA_warn(__name__)


def grab_fpy(build_dir="", file_out="wimsd-fpyield.html"):
    """Grabs the WIMS fission product yields from the IAEA website"""
    build_filename = os.path.join(build_dir, file_out)
    local_filename = os.path.join(os.path.dirname(__file__), file_out)

    if os.path.exists(local_filename):
        shutil.copy(local_filename, build_filename)
        return

    nist = urllib2.urlopen("https://www-nds.iaea.org/wimsd/fpyield.htm")
    with open(build_filename, "w") as f:
        f.write(nist.read())


class Parser(HTMLParser):
    """Parser for WIMSD fission product yield files."""

    def __init__(self, *args, **kwargs):
        # HTMLParser is not a new style class, no super()
        HTMLParser.__init__(self, *args, **kwargs)
        self._currrow = []
        self._fromnucs = None
        self.fission_product_yields = []

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self._currrow = []

    def handle_endtag(self, tag):
        if tag == "tr":
            row = self._currrow
            self._currrow = []
            if len(row) == 7 and row[0] == "Fission" and row[1] == "product":
                self._fromnucs = [nucname.id(x) for x in row[-4:]]
                return
            if len(row) != 6:
                return
            if row[0].endswith("FP"):
                return
            tonuc = nucname.id(row[0].split("-", 1)[1])
            self.fission_product_yields += zip(
                self._fromnucs, [tonuc] * 4, map(float, row[-4:])
            )

    def handle_data(self, data):
        data = "".join(data.split())
        if len(data) == 0:
            return
        self._currrow.append(data.strip())


fpy_dtype = np.dtype(
    [
        ("from_nuc", "i4"),
        ("to_nuc", "i4"),
        ("yields", "f8"),
    ]
)


def parse_fpy(build_dir):
    """Converts fission product yeild data to a numpy array."""
    build_filename = os.path.join(build_dir, "wimsd-fpyield.html")
    with open(build_filename, "r") as f:
        raw_data = f.read()
    parser = Parser()
    parser.feed(raw_data)
    data = sorted(parser.fission_product_yields)
    data = np.array(data, dtype=fpy_dtype)
    return data


def make_fpy_table(nuc_data, build_dir=""):
    """Adds the neutron scattering lengths to the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory the html files in.
    """
    yields = parse_fpy(build_dir)
    db = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)
    if not hasattr(db.root, "neutron"):
        neutron_group = db.create_group("/", "neutron", "Neutron Data")
    fpy_table = db.create_table(
        "/neutron/",
        "wimsd_fission_products",
        yields,
        "WIMSD Fission Product Yields, fractions [unitless]",
    )
    fpy_table.flush()
    db.close()


def make_fpy(args):
    """Controller function for WIMS fission products."""
    nuc_data, build_dir = args.nuc_data, args.build_dir

    # Check that the table exists
    with tb.open_file(nuc_data, "a", filters=BASIC_FILTERS) as f:
        if hasattr(f.root, "neutron") and hasattr(
            f.root.neutron, "wimsd_fission_products"
        ):
            print(
                "skipping WIMSD fission product yield table creation; "
                "already exists."
            )
            return

    # Grab the raw data
    print("Grabbing WIMSD fission product yield data.")
    grab_fpy(build_dir)

    # Make scatering table once we have the data
    print("Making WIMSD fission product yield table.")
    make_fpy_table(nuc_data, build_dir)
