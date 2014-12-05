"""Module handles the construction of a reference materials library in nuc_data.h5.
This currently consists to natural element materials and those coming from PNNL's
`Materials Compendium`_.

.. _Materials Compendium: http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev1.pdf
"""

from __future__ import print_function
import re
import os
import csv
from itertools import takewhile, groupby
from warnings import warn
from pyne.utils import QAWarning

import tables as tb
import numpy as np

from pyne import nucname
from pyne.data import natural_abund, natural_abund_map
from pyne.material import Material

warn(__name__ + " is not yet QA compliant.", QAWarning)

nucids = set()
natural_abund("H1")  # initialize natural_abund_map


# Make a dictionary that represents elements as dicts of their isotopes
def make_elements():
    """Make natural elemental materials based on isotopic abundances.

    Returns
    -------
    elts : [pyne.material.Material]
        Natural elements as materials.
    """
    sorted_abunds = sorted(list(natural_abund_map.items()))
    grouped_abunds = groupby(sorted_abunds, lambda x: nucname.zzzaaa(x[0])//1000)
    elts = (Material(dict(abunds), metadata={"name": nucname.name(zz)})
            for zz, abunds in grouped_abunds)
    return elts


# Parses data from .csv
def grab_materials_compendium(location='materials_compendium.csv'):
    """Parses data from a materials compendium csv file.

    Parameters
    ----------
    location : str
        The file to read in compendium from.

    Returns
    -------
    mats : [pyne.material.Material]
        The materials in the compendium.
    """
    with open(location, 'r', newline='') as f:
        lines = csv.reader(f, delimiter=',', quotechar='"')
        lines = list(filter(is_comp_matname_or_density, lines))
        mats = parse_materials([], lines)
        return mats


def is_comp_matname_or_density(line):
    """Detect composition, material name, or density lines.

    Parameters
    ----------
    line : [str]
        The input line.

    Returns
    -------
    result : bool
        True if the input line has composition, material name, or density data.
        False otherwise.
    """
    if not line[0]:
        return False
    if line[0] == "Density (g/cm3) =":
        return True
    if re.match(r'\d+. +$|[A-Za-z]{1,2}-?(\d{1,3})?$', line[0]):
        return True
    return False


def elem_line_to_mat(line):
    """Take an element line and turn it into a PyNE material

    Parameters
    ==========
    line : [str]
        The input line.

    Returns
    =======
    mat : pyne.material.Material
        The material represented by this line, weighted by mass fraction.
    """
    name = line[0]
    mass_frac = float(line[3])
    try:
        za = int(line[1])
        # maybe we should check if za % 1000 == 0?
        mat = Material({nucname.id(za): mass_frac})
    except ValueError:
        za = None
        mat = expand_elt_to_mat(name) * mass_frac
    return mat


def expand_elt_to_mat(elt_id):
    """Expands an element into a material based on its natural abundances.

    Parameters
    ----------
    elt_id : int or str
        The element you wish to expand.

    Returns
    -------
    mat : pyne.material.Material
        A material with the isotopic abundances.
    """
    if nucname.anum(elt_id) != 0:
        raise ValueError("Expected an element, got a specific nuclide instead.")
    elt_z = nucname.zzzaaa(elt_id) // 1000
    nucs = ((nuc, abund) for (nuc, abund) in natural_abund_map.items()
            if nucname.zzzaaa(nuc) // 1000 == elt_z and nucname.anum(nuc) != 0)
    return Material(dict(nucs))


def parse_materials(mats, lines):
    """Take first material from ``lines`` and append to ``mats``.

    Parameters
    ----------
    mats : [pyne.material.Material]
        The growing list of materials.
    lines: [[str]]
        The shrinking list of lines.
    """
    if len(lines) == 0:
        return mats
    material_lines = list(takewhile(lambda l: re.match(r"^\d+. +", l[0]) is None, lines[2:]))
    material_length = len(material_lines) + 2
    mat = sum((elem_line_to_mat(l) for l in material_lines))
    mat.density = float(lines[1][2])
    mat.metadata = {"name": lines[0][1]}
    mat.normalize()
    mats.append(mat)
    return parse_materials(mats, lines[material_length:])

# Writes to file
def make_materials_compendium(nuc_data, mats, elts):
    """Adds materials compendium to nuc_data.h5."""
    # open nuc_data, make nuc_zz an array
    filters = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)
    with tb.openFile(nuc_data, 'r+', filters=filters) as f:
        f.createGroup('/', 'material_library')
        f.createArray('/material_library', 'nucid', np.array(sorted(nucids)))

    for elt in elts:
        elt.write_hdf5(nuc_data, datapath="/material_library/materials",
                           nucpath="/material_library/nucid", chunksize=70)

    # Writes materials from mats to file, and names them.
    for mat in mats:
        mat.write_hdf5(nuc_data, datapath="/material_library/materials",
                       nucpath="/material_library/nucid", chunksize=70)


def make_materials_library(args):
    """Controller function for adding materials library."""
    nuc_data = args.nuc_data
    if os.path.exists(nuc_data):
        with tb.openFile(nuc_data, 'r') as f:
            if '/material_library' in f:
                print("skipping materials library data table creation; already exists.")
                return

    print("Making the elements...")
    elts = make_elements()
    print("Grabbing materials compendium...")
    mats = grab_materials_compendium(os.path.join(os.path.split(__file__)[0],
                                                  'materials_compendium.csv'))
    print("Making materials library...")
    make_materials_compendium(nuc_data, mats, elts)

if __name__ == "__main__":
    make_elements()
        
