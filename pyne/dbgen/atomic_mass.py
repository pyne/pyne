"""This module provides a way to grab and store raw data for atomic mass."""
from __future__ import print_function
import os
import re
import pkgutil
from pyne.utils import QA_warn

import numpy as np
import tables as tb

from pyne import nucname
from pyne.dbgen.api import BASIC_FILTERS
from pyne.dbgen.isotopic_abundance import get_isotopic_abundances

QA_warn(__name__)

# Note that since ground state and meta-stable isotopes are of the same atomic mass,
# the meta-stables have been discluded from the following data sets.

MASS_FILE = "mass.mas16"


def copy_atomic_mass_adjustment(build_dir=""):
    """Copies the atomic mass evaluation originally from the Atomic Mass Data
    Center.  These are courtesy of Georges Audi and Wang Meng via a private
    communication, November 2012."""

    if os.path.exists(os.path.join(build_dir, MASS_FILE)):
        return

    mass = pkgutil.get_data("pyne.dbgen", MASS_FILE)
    with open(os.path.join(build_dir, MASS_FILE), "wb") as f:
        f.write(mass)


# Note, this regex specifically leaves our free neutrons
# amdc_regex = re.compile('[ \d-]*? (\d{1,3})[ ]{1,4}(\d{1,3}) [A-Z][a-z]? .*? (\d{1,3}) ([ #.\d]{10,11}) ([ #.\d]{1,10})[ ]*?$')
amdc_regex = re.compile(
    "[ \d-]*? (\d{1,3})[ ]{1,4}(\d{1,3}) [A-Z][a-z]? .*? (\d{1,3}) ([ #.\d]{5,12}) ([ #.\d]+)[ ]*?$"
)


def parse_atomic_mass_adjustment(build_dir=""):
    """Parses the atomic mass adjustment data into a list of tuples of
    the nuclide, atomic mass, and error."""
    f = open(os.path.join(build_dir, MASS_FILE), "r")

    atomic_masses = []

    for line in f:
        m = amdc_regex.search(line)
        if m is None:
            continue

        nuc = (10000000 * int(m.group(1))) + (10000 * int(m.group(2)))
        mass = float(m.group(3)) + 1e-6 * float(m.group(4).strip().replace("#", ""))
        error = 1e-6 * float(m.group(5).strip().replace("#", ""))

        atomic_masses.append((nuc, mass, error))

    f.close()

    return atomic_masses


atomic_mass_desc = {
    "nuc": tb.IntCol(pos=1),
    "mass": tb.FloatCol(pos=2),
    "error": tb.FloatCol(pos=3),
    "abund": tb.FloatCol(pos=4),
}

atomic_mass_dtype = np.dtype(
    [
        ("nuc", int),
        ("mass", float),
        ("error", float),
        ("abund", float),
    ]
)


def make_atomic_mass_table(nuc_data, build_dir=""):
    """Makes an atomic mass table in the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place html files in.
    """
    # Grab raw data
    atomic_abund = get_isotopic_abundances()
    atomic_masses = parse_atomic_mass_adjustment(build_dir)

    A = {}

    # Add normal isotopes to A
    for nuc, mass, error in atomic_masses:
        if nuc in atomic_abund:
            A[nuc] = nuc, mass, error, atomic_abund[nuc]
        else:
            A[nuc] = nuc, mass, error, 0.0

    # Add naturally occuring elements
    for element in nucname.name_zz:
        nuc = nucname.id(element)
        A[nuc] = nuc, 0.0, 0.0, 0.0

    for nuc, abund in atomic_abund.items():
        zz = nucname.znum(nuc)
        element_zz = nucname.id(zz)
        element = nucname.zz_name[zz]

        _nuc, nuc_mass, _error, _abund = A[nuc]
        elem_zz, elem_mass, _error, _abund = A[element_zz]

        new_elem_mass = elem_mass + (nuc_mass * abund)
        A[element_zz] = element_zz, new_elem_mass, 0.0, float(0.0 < new_elem_mass)

    A = sorted(A.values(), key=lambda x: x[0])

    # Open the HDF5 File
    kdb = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Make a new the table
    Atable = kdb.create_table(
        "/",
        "atomic_mass",
        atomic_mass_desc,
        "Atomic Mass Data [amu]",
        expectedrows=len(A),
    )
    Atable.append(A)

    # Ensure that data was written to table
    Atable.flush()

    # Close the hdf5 file
    kdb.close()


def make_atomic_mass(args):
    """Controller function for adding atomic_mass."""
    nuc_data, build_dir = args.nuc_data, args.build_dir

    if os.path.exists(nuc_data):
        with tb.open_file(nuc_data, "r") as f:
            if hasattr(f.root, "atomic_mass"):
                print("skipping atomic mass data table creation; already exists.")
                return

    # Then grab mass data
    print("Copying AME 2016 atomic mass data.")
    copy_atomic_mass_adjustment(build_dir)

    # Make atomic mass table once we have the array
    print("Making atomic mass data table.")
    make_atomic_mass_table(nuc_data, build_dir)
