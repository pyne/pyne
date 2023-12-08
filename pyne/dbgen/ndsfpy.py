"""This module provides a way to grab and store raw data for fission product yeilds
from the NDS library at the IAEA. For more information, please visit their website:
https://www-nds.iaea.org/sgnucdat/index.htm or
https://www-nds.iaea.org/sgnucdat/c2.htm. Please contact the NDS at
online@iaeand.iaea.org with questions about the data itself.

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
from __future__ import print_function, division
import os
import shutil
from pyne.utils import QA_warn

try:
    import urllib.request as urllib2
except ImportError:
    import urllib2

import numpy as np
import numpy.lib.recfunctions
import tables as tb

from pyne import nucname
from pyne.dbgen.api import BASIC_FILTERS

QA_warn(__name__)


def readtable(i, spdat):
    """
    Reads in a set of 5 html tables and returns corresponding yield data
    """
    parent = getdata(i, spdat)[0]
    pfinal = (parent.split("<strong>")[1]).split("</strong>")[0]
    pid = conv_to_id(pfinal)
    fpdata = getdata(i + 1, spdat)
    dt = np.dtype(
        [
            ("from_nuc", "i4"),
            ("to_nuc", "i4"),
            ("yield_thermal", float),
            ("yield_thermal_err", float),
            ("yield_fast", float),
            ("yield_fast_err", float),
            ("yield_14MeV", float),
            ("yield_14MeV_err", float),
        ]
    )
    dfinal = np.zeros((len(fpdata),), dtype=dt)
    for index, item in enumerate(fpdata):
        dfinal[index]["from_nuc"] = pid
        dfinal[index]["to_nuc"] = conv_to_id(item)
    thermaldata = getdata(i + 2, spdat)
    for index, item in enumerate(thermaldata):
        dat, err = conv_to_num(item)
        dfinal[index]["yield_thermal"] = dat
        dfinal[index]["yield_thermal_err"] = err
    fastdata = getdata(i + 3, spdat)
    for index, item in enumerate(fastdata):
        dat, err = conv_to_num(item)
        dfinal[index]["yield_fast"] = dat
        dfinal[index]["yield_fast_err"] = err
    dtdata = getdata(i + 4, spdat)
    for index, item in enumerate(dtdata):
        dat, err = conv_to_num(item)
        dfinal[index]["yield_14MeV"] = dat
        dfinal[index]["yield_14MeV_err"] = err
    return dfinal


def conv_to_id(nuc):
    """Converts html nuclide names to nuclide ids"""
    parts = nuc.split("-")
    return nucname.id(parts[1] + parts[2])


def conv_to_num(dstring):
    """Converts html number and error to floats"""
    if dstring == "-":
        return 0, 0
    dat, err = dstring.split("&plusmn;")
    if "<sup>" in dat:
        dat = parse_num(dat)
    else:
        dat = float(dat)
    if "<sup>" in err:
        err = parse_num(err)
    else:
        err = float(err)
    return dat, err


def parse_num(dst):
    """Converts html numbers with exponents to floats"""
    nums = dst.split("x")
    base = float(nums[0])
    exp = (nums[1].split("<sup>")[1]).split("</sup>")[0]
    return base * 10 ** float(exp)


def getpoint(line):
    """Gets data entries from html lines"""
    spline = line.split('<tr><td class="xl28b">&nbsp;&nbsp;')
    if len(spline) > 1:
        data = spline[1].split("</td></tr>")[0]
    else:
        data = None
    return data


def getdata(i, spdat):
    """Gets the data from the nds html table"""
    lines = spdat[i].splitlines()
    dlist = []
    for line in lines:
        d = getpoint(line)
        if d is None:
            continue
        dlist.append(d)
    return dlist


def make_fpy_table(nuc_data, build_dir=""):
    """Adds the NDS fission yields to the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    """
    build_filename = os.path.join(build_dir, "nds-fpyield.html")
    with open(build_filename, "rb") as f:
        raw_data = f.read().decode("iso-8859-1")
    spdat = raw_data.split("<table>")
    alldata = []
    for i in range(1, 31, 5):
        alldata.append(readtable(i, spdat))
    alldata = numpy.lib.recfunctions.stack_arrays(alldata, asrecarray=True)
    db = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)
    if not hasattr(db.root, "neutron"):
        neutron_group = db.create_group("/", "neutron", "Neutron Data")
    fpy_table = db.create_table(
        "/neutron/",
        "nds_fission_products",
        alldata,
        "NDS Fission Product Yields, percent [unitless]",
    )
    fpy_table.flush()
    db.close()


def grab_fpy(build_dir="", file_out="nds-fpyield.html"):
    """Grabs the NDS fission product yields from the IAEA website"""
    build_filename = os.path.join(build_dir, file_out)
    local_filename = os.path.join(os.path.dirname(__file__), file_out)

    if os.path.exists(local_filename):
        shutil.copy(local_filename, build_filename)
        return

    nist = urllib2.urlopen("https://www-nds.iaea.org/sgnucdat/c3.htm")
    with open(build_filename, "wb") as f:
        f.write(nist.read())


def make_fpy(args):
    """Controller function for NDS fission products."""
    nuc_data, build_dir = args.nuc_data, args.build_dir
    # Check that the table exists
    with tb.open_file(nuc_data, "a", filters=BASIC_FILTERS) as f:
        if hasattr(f.root, "neutron") and hasattr(
            f.root.neutron, "nds_fission_products"
        ):
            print(
                "skipping NDS fission product yield table creation; " "already exists."
            )
            return
    print("Grabbing NDS fission product yield data.")
    grab_fpy(build_dir)

    print("Making NDS fission product yield table.")
    make_fpy_table(nuc_data, build_dir)
