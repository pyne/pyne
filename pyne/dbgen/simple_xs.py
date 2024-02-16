"""This module provides a way to grab and store simple cross sections from KAERI."""
from __future__ import print_function
import os
from pyne.utils import QA_warn

try:
    import urllib.request as urllib
except ImportError:
    import urllib
from zipfile import ZipFile

import numpy as np
import tables as tb

from .. import nucname
from ..utils import to_barns
from .api import BASIC_FILTERS
from .kaeri import grab_kaeri_nuclide, parse_for_all_isotopes

QA_warn(__name__)


def grab_kaeri_simple_xs(build_dir=""):
    """Grabs the KAERI files needed for the simple cross sections table,
    if not already present.

    Parameters
    ----------
    build_dir : str
        Major directory to place html files in. 'KAERI/' will be appended.
    """
    zip_path = os.path.join(build_dir, "kaeri.zip")
    zip_url = "http://data.pyne.io/kaeri.zip"
    if not os.path.exists(zip_path):
        print("  grabbing {0} and placing it in {1}".format(zip_url, zip_path))
        urllib.urlretrieve(zip_url, zip_path)
        try:
            zf = ZipFile(zip_path)
            for name in zf.namelist():
                if not os.path.exists(os.path.join(build_dir, name)):
                    print("    extracting {0} from {1}".format(name, zip_path))
                    zf.extract(name, build_dir)
        finally:
            zf.close()

    # Add kaeri to build_dir
    build_dir = os.path.join(build_dir, "KAERI")
    try:
        os.makedirs(build_dir)
    except OSError:
        pass
    already_grabbed = set(os.listdir(build_dir))

    # Grab and parse elemental summary files.
    all_nuclides = set()
    for element in nucname.name_zz.keys():
        htmlfile = element.upper() + ".html"
        if htmlfile not in already_grabbed:
            grab_kaeri_nuclide(element.upper(), build_dir)
        all_nuclides = all_nuclides | parse_for_all_isotopes(
            os.path.join(build_dir, htmlfile)
        )

    # Grab nuclide XS summary files
    for nuc in sorted(all_nuclides):
        nuc = nucname.name(nuc)
        htmlfile = nuc.upper() + "_2.html"
        if htmlfile not in already_grabbed:
            grab_kaeri_nuclide(nuc.upper(), build_dir, 2)


simple_xs_channels = {
    "sigma_t": "Total Cross Section",
    "sigma_e": "Elastic Scattering Cross Section",
    "sigma_i": "Total Inelastic Cross Section",
    "sigma_2n": "(n,2n) Cross Section",
    "sigma_3n": "(n,3n) Cross Section",
    "sigma_4n": "(n,4n) Cross Section",
    "sigma_f": "Total Fission Cross Section",
    "sigma_gamma": "Radiative Capture Cross Section",
    "sigma_alpha": "(n,alpha) Cross Section",
    "sigma_proton": "(n,p) Cross Section",
    "sigma_deut": "(n,d) Cross Section",
    "sigma_trit": "(n,t) Cross Section",
}


simple_xs_energy = {
    "thermal": "at 0.0253 eV",
    "thermal_maxwell_ave": "Maxwell avg. at 0.0253 eV",
    "resonance_integral": "Resonance integral",
    "fourteen_MeV": "at 14 MeV",
    "fission_spectrum_ave": "Fission spectrum avg.",
}


simple_xs_dtype = np.dtype(
    [
        ("nuc", int),
        ("sigma_t", float),
        ("sigma_s", float),
        ("sigma_e", float),
        ("sigma_i", float),
        ("sigma_a", float),
        ("sigma_gamma", float),
        ("sigma_f", float),
        ("sigma_alpha", float),
        ("sigma_proton", float),
        ("sigma_deut", float),
        ("sigma_trit", float),
        ("sigma_2n", float),
        ("sigma_3n", float),
        ("sigma_4n", float),
    ]
)


def get_xs_from_file(filename, eng, chan):
    """Parses out a cross section from a KAERI file.

    Parameters
    ----------
    filename : str
        Local path to a KAERI neutron cross section summary html file.
    eng : str
        Energy flag to find this cross section for.  (Must be key
        of simple_xs_energy dictionary).
    chan : str
        Cross section (interaction channel) to find.  (Must be key
        of simple_xs_channels dict).

    Returns
    -------
    data : float
        Microscopic cross section in [barns].
    """
    with open(filename, "r") as f:
        in_channel = False
        for line in f:
            if simple_xs_channels[chan] in line:
                in_channel = True

            if in_channel and ("<li>" + simple_xs_energy[eng] in line):
                du = line.partition("=")[2].split()
                data = float(du.pop(0))
                unit = ""
                for u in du:
                    unit = unit + u
                unit = unit.partition("\\")[0]
                data = to_barns(data, unit)
                return data

            elif in_channel and ("</ul>" in line):
                # XS not defined for this energy, returning zero
                return 0.0

    # If the specific XS was not found in this file, return zero
    return 0.0


def parse_simple_xs(build_dir=""):
    """Builds and returns a dictionary from cross-section types to nuclides."""
    build_dir = os.path.join(build_dir, "KAERI")

    # Grab and parse elemental summary files.
    all_nuclides = set()
    for element in nucname.name_zz.keys():
        htmlfile = element.upper() + ".html"
        all_nuclides = all_nuclides | parse_for_all_isotopes(
            os.path.join(build_dir, htmlfile)
        )
    all_nuclides = sorted([nucname.id(nuc) for nuc in all_nuclides])
    energy_tables = dict(
        [
            (eng, np.zeros(len(all_nuclides), dtype=simple_xs_dtype))
            for eng in simple_xs_energy.keys()
        ]
    )

    # Loop through species
    for i, nuc in enumerate(all_nuclides):
        nuc_name = nucname.name(nuc).upper()
        filename = os.path.join(build_dir, nuc_name + "_2.html")

        # Loop through all energy types
        for eng in simple_xs_energy:
            energy_tables[eng]["nuc"][i] = nuc

            # Loop trhough reactions
            for chan in simple_xs_channels:
                energy_tables[eng][chan][i] = get_xs_from_file(filename, eng, chan)

    for eng in simple_xs_energy:
        # Store only non-trivial entries
        channels_list = list(simple_xs_channels.keys())
        mask = (
            energy_tables[eng][channels_list]
            != np.zeros(1, dtype=simple_xs_dtype)[channels_list]
        )
        energy_tables[eng] = energy_tables[eng][mask]

        # Calculate some xs
        energy_tables[eng]["sigma_s"] = (
            energy_tables[eng]["sigma_e"] + energy_tables[eng]["sigma_i"]
        )

        energy_tables[eng]["sigma_a"] = (
            energy_tables[eng]["sigma_gamma"]
            + energy_tables[eng]["sigma_f"]
            + energy_tables[eng]["sigma_alpha"]
            + energy_tables[eng]["sigma_proton"]
            + energy_tables[eng]["sigma_deut"]
            + energy_tables[eng]["sigma_trit"]
            + energy_tables[eng]["sigma_2n"]
            + energy_tables[eng]["sigma_3n"]
            + energy_tables[eng]["sigma_4n"]
        )
    return energy_tables


def make_simple_xs_tables(nuc_data, build_dir=""):
    """Make the simple cross section tables.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place html files in.
    """
    # Grab raw data
    simple_xs_tables = parse_simple_xs(build_dir)

    # Open the HDF5 File
    db = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Create neutron group
    if not hasattr(db.root, "neutron"):
        neutron_group = db.create_group("/", "neutron", "Neutron Interaction Data")

    # Create simple_xs Group
    if not hasattr(db.root.neutron, "simple_xs"):
        simple_xs_group = db.create_group(
            "/neutron", "simple_xs", "Simple Neutron Cross Section Data"
        )

    # Create tables for every energy
    for eng, eng_flag in simple_xs_energy.items():
        simple_xs_table = db.create_table(
            simple_xs_group,
            eng,
            np.empty(0, dtype=simple_xs_dtype),
            "{0} [barns]".format(eng_flag.capitalize()),
            expectedrows=len(simple_xs_tables[eng]),
        )
        simple_xs_table.append(simple_xs_tables[eng])
        simple_xs_table.flush()

    # Close the hdf5 file
    db.close()


def make_simple_xs(args):
    """Controller function for adding basic cross section data."""
    nuc_data, build_dir = args.nuc_data, args.build_dir

    with tb.open_file(nuc_data, "a", filters=BASIC_FILTERS) as f:
        if hasattr(f.root, "neutron") and hasattr(f.root.neutron, "simple_xs"):
            print("skipping simple XS data table creation; already exists.")
            return

    # First grab the atomic abundance data
    print("Grabbing neutron summary files from KAERI")
    grab_kaeri_simple_xs(build_dir)

    # Make simple table once we have the array
    print("Making simple cross section data tables")
    make_simple_xs_tables(nuc_data, build_dir)
