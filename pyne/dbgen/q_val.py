"""
Module allows the grabbing of q_values (energy per disintegration) for the
calculation of decay heat. This currently consists of the nuclide, it's
q_value, and the percent of q coming from gammas. This data is from
'ORIGEN-S DECAY DATA LIBRARY AND HALF-LIFE UNCERTAINTIES'
(http://web.ornl.gov/~webworks/cppr/y2001/rpt/97914.pdf)
"""
from __future__ import print_function
import csv
import os
from pyne.utils import QA_warn

import numpy as np
import tables as tb

from pyne import nucname
from pyne.api import nuc_data
from pyne.dbgen.api import BASIC_FILTERS

QA_warn(__name__)

# Parses data from .csv
def grab_q_values(fname):
    """Parses data from three q_val csv files.

    Parameters
    ----------
    fname : str
        Name of q_value file.
    """

    # Create list
    all_q_values = []

    # Open .csv files and parses them
    with open(os.path.join(os.path.dirname(__file__), fname), "r") as f:
        reader = csv.reader(f)
        for row in reader:
            entry = read_row(row)
            all_q_values.append(entry)

    return all_q_values


# Grabs data row by row
def read_row(row):
    """Returns a list of the format [int, float, float] for each nuclide.

    Parameters
    ----------
    row : tuple
        One entry in a q_val file.
    """

    # Create tuple
    entry = ()

    # Evaluate each component of the given row
    if row[0] == "Nuclide" or len(row[0].strip()) == 0:
        return
    nuclide = nucname.id("".join(row[0:2]).replace(" ", ""))
    if len(row[2]) == 0:
        q_val = 0.0
    else:
        q_val = float(row[2])
    if len(row[3]) == 0:
        gamma_frac = 0.0
    else:
        gamma_frac = float(row[3])
    entry = (nuclide, q_val, gamma_frac)

    return entry


# Sorts and filters list of q_values
def format_q_values(all_q_values):
    """Filters the q_value data for multiple entries then sorts the nuclides.

    Parameters
    ----------
    all_q_values : list of tuples
        Array of q_values for all nuclides.
    """

    # Create list
    d_all_q_values = []
    distinct_all_q_values = []

    # Ensure only one entry per nuclide then sort in order of nuclide
    for nuclide in all_q_values:
        if nuclide is not None:
            d_all_q_values.append(nuclide)
    distinct_all_q_values = list(set(d_all_q_values))
    distinct_all_q_values.sort(key=lambda nucid: nucid[0])

    return distinct_all_q_values


# Write q_value table to file
def make_q_value_table(all_q_values, nuc_data, build_dir=""):
    """Adds q_value table to the nuc_data.h5 library.

    Parameters
    ----------
    all_q_values: list of tuples
        Array of q_values for all nuclides.
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place q_value files in.
    """

    # Sort and filter the q_values, make into list of tuples
    distinct_all_q_values = format_q_values(all_q_values)

    # Define data type
    qv_dtype = np.dtype(
        [
            ("nuc", np.int32),
            ("q_val", np.float64),
            ("gamma_frac", np.float64),
        ]
    )

    # Convert to numpy array
    q_value_array = np.array(distinct_all_q_values, dtype=qv_dtype)

    # Open the hdf5 file
    nuc_file = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Make the group if it's not there
    if not hasattr(nuc_file.root, "decay"):
        nuc_file.create_group("/", "decay", "ENSDF Decay data")

    # Make a new table
    q_value_table = nuc_file.create_table(
        "/decay",
        "q_values",
        q_value_array,
        "Nuclide, Q_value [MeV per disintegration], Fraction of Q that comes from gammas",
    )

    # Ensure that data was written to table
    q_value_table.flush()

    # Close the hdf5 file
    nuc_file.close()


def make_q_value(args):
    """Controller function for adding q-values"""
    nuc_data, build_dir = args.nuc_data, args.build_dir
    if os.path.exists(nuc_data):
        with tb.open_file(nuc_data, "r") as f:
            if "/decay/q_values" in f:
                print("skipping q_value table creation; already exists.")
                return

    # Grab the q_values
    print("Grabbing q_values...")
    q_value_files = [
        "q_val_actinides.csv",
        "q_val_fissionproducts.csv",
        "q_val_light.csv",
    ]
    all_q_values = []
    for fname in q_value_files:
        all_q_values += grab_q_values(fname)

    # Make the q_value table and write to file
    print("Making q_value table...")
    make_q_value_table(all_q_values, nuc_data, build_dir)
