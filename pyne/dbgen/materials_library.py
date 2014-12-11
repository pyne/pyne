"""Module handles the construction of a reference materials library in nuc_data.h5.
This currently consists to natural element materials and those coming from PNNL's
`Materials Compendium`_.

.. _Materials Compendium: http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev1.pdf
"""

from __future__ import print_function
import re
import os
import csv
import sys
from itertools import takewhile, groupby
from warnings import warn
from pyne.utils import QAWarning

import tables as tb

from pyne import nucname
from pyne.data import natural_abund, natural_abund_map
from pyne.material import Material, MaterialLibrary

warn(__name__ + " is not yet QA compliant.", QAWarning)


# Make a dictionary that represents elements as dicts of their isotopes
def make_elements():
    """Make natural elemental materials based on isotopic abundances.

    Returns
    -------
    eltsdict : {str: pyne.material.Material}
        Natural elements as materials.
    """
    natural_abund("H1")  # initialize natural_abund_map
    sorted_abunds = sorted(list(natural_abund_map.items()))
    grouped_abunds = groupby(sorted_abunds, lambda x: nucname.zzzaaa(x[0])//1000)
    # filter out 111, 113, 115, 117, 118 - the ones with no names
    elts = (Material(dict(abunds), metadata={"name": nucname.name(zz)})
            for zz, abunds in grouped_abunds if zz in nucname.zz_name.keys())
    eltsdict = dict(((elt.metadata["name"], elt) for elt in elts))
    return eltsdict


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
    natural_abund("H1")  # initialize natural_abund_map
    try:
        f = open(location, 'r', newline='')
    except TypeError:  # python2
        f = open(location, 'rb')
    reader = csv.reader(f, delimiter=',', quotechar='"')
    lines = list(filter(is_comp_matname_or_density, reader))
    mats = parse_materials({}, lines)
    f.close()
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


def parse_materials(mats, lines):
    """Take first material from ``lines`` and append to ``mats``.

    Parameters
    ----------
    mats : {str: pyne.material.Material}
        The growing dict of materials.
    lines: [[str]]
        The shrinking list of lines.
    """
    if len(lines) == 0:
        return mats
    material_lines = list(takewhile(lambda l: re.match(r"^\d+. +", l[0]) is None, lines[2:]))
    material_length = len(material_lines) + 2
    # mat = sum((elem_line_to_mat(l) for l in material_lines))
    mat = sum((Material({l[0]: float(l[3])}) for l in material_lines))
    mat.density = float(lines[1][2])
    name = lines[0][1]
    mat.metadata = {"name": name}
    mat.normalize()
    mat = mat.expand_elements()
    mats.update({name: mat})
    return parse_materials(mats, lines[material_length:])


# Writes to file
def make_materials_compendium(nuc_data, matslib):
    """Adds materials compendium to nuc_data.h5."""
    matslib.write_hdf5(nuc_data, datapath="/material_library/materials",
                       nucpath="/material_library/nucid")


def make_materials_library(args):
    """Controller function for adding materials library."""
    nuc_data = args.nuc_data
    if os.path.exists(nuc_data):
        with tb.openFile(nuc_data, 'r') as f:
            if '/material_library' in f:
                print("skipping materials library data table creation; already exists.")
                return

    print("Making the elements...")
    eltsdict = make_elements()
    matslib = MaterialLibrary(eltsdict)
    print("Grabbing materials compendium...")
    matsdict = grab_materials_compendium(os.path.join(os.path.split(__file__)[0],
                                                      'materials_compendium.csv'))
    matslib.update(matsdict)
    print("Making materials library...")
    make_materials_compendium(nuc_data, matslib)
