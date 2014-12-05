"""Module handles the construction of a reference materials library in nuc_data.h5.
This currently consists to natural element materials and those coming from PNNL's
`Materials Compendium`_.

.. _Materials Compendium: http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev1.pdf
"""

from __future__ import print_function
import re
import os
import csv
from itertools import takewhile, islice, filterfalse
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
    """Make natural elemental materials based on isotopic abundances."""
    elemental_mats = {}
    for name, zz in nucname.name_zz.items():
        elemental_mats[name] = {}
    for nuc, abund in natural_abund_map.items():
        nucid = nucname.id(nuc)
        anum = nucname.anum(nucid)
        if 0 == anum or abund == 0.0:
            continue
        znum = nucname.znum(nuc)
        if znum not in nucname.zz_name:
            continue
        name = nucname.zz_name[znum]
        elemental_mats[name][nucid] = abund
        nucids.add(nucid)
    return elemental_mats


# Parses data from .csv
def grab_materials_compendium(location='materials_compendium.csv'):
    """Parses data from a materials compendium csv file."""
    with open(location, 'r', newline='') as f:
        lines = csv.reader(f, delimiter=',', quotechar='"')
        lines = list(filter(is_comp_matname_or_density, lines))
        mats = parse_materials([], lines)
        return mats


def is_comp_matname_or_density(line):
    return ((line[0] == "Density (g/cm3) =" or
             re.match(r'\d+. +$|[A-Za-z]{1,2}-?(\d{1,3})?$', line[0])))


def elem_line_to_mat(line):
    """Take an element line and turn it into a PyNE material

    Parameters
    ==========
    line : [str]
        The input line.

    Returns
    =======
    mat : pyne.Material
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
    try:
        name = lines[0][1]
    except IndexError:
        return mats
    density = float(lines[1][2])
    material_lines = list(takewhile(lambda l: re.match(r"^\d+. +", l[0]) is None, lines[2:]))
    material_length = len(material_lines) + 2
    mat = sum((elem_line_to_mat(l) for l in material_lines))
    mat.density = density
    mat.metadata = {"name": name}
    mat.normalize()
    mats.append(mat)
    # return parse_materials(mats, islice(lines, material_length))
    return parse_materials(mats, lines[material_length:])

# Writes to file
def make_materials_compendium(nuc_data, mats, elts):
    """Adds materials compendium to nuc_data.h5."""
    # open nuc_data, make nuc_zz an array
    filters = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)
    with tb.openFile(nuc_data, 'r+', filters=filters) as f:
        f.createGroup('/', 'material_library')
        f.createArray('/material_library', 'nucid', np.array(sorted(nucids)))

    # Writes elements for which we have compositional data to file
    for zz in elts:
        if 0 == len(elts[zz]):
            continue
        element = Material(elts[zz], mass=1.0,
                           metadata={'name': nucname.name(zz)})
        element.write_hdf5(nuc_data, datapath="/material_library/materials",
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

    # First make the elements
    print("Making the elements...")
    elts = make_elements()

    # Then grab the materials compendium
    print("Grabbing materials compendium...")
    mats = grab_materials_compendium(os.path.join(os.path.split(__file__)[0],
                                                  'materials_compendium.csv'))

    # Make atomic mass table once we have the array
    print("Making materials library...")
    make_materials_compendium(nuc_data, mats, elts)


if __name__ == "__main__":
    elts = make_elements()
    mats = grab_materials_compendium()
    import pdb; pdb.set_trace()
        
