"""Module handles the construction of a reference materials library in nuc_data.h5.
This currently consists to natural element materials and those coming from PNNL's
`Materials Compendium`_.

.. _Materials Compendium: http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev1.pdf
"""

import csv    
import re
from pyne.material import Material
import tables as tb
import numpy as np
from pyne import nucname
import os
from pyne.data import natural_abund, natural_abund_map

elemental_mats = {}
names = []
nucids = set()
mats = []
densities = []  

# Make a dictionary that represents elements as dicts of their isotopes
def make_elements():
    """Makes natural elemental materials based on isotopic abundances."""
    habund = natural_abund('H')
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
    
# Parses data from .csv
def grab_materials_compendium(location = 'materials_compendium.csv'):
    """Parses data from a materials compendium csv file."""
    # grabs name from starting row, starts a new dictionary for composition
    def starting_row(row):
        if re.match('\d{1,3}\.  ', row[0]):
            #print row
            name = row[1]
            names.append(name)
            composition.clear()
        else:
            pass

    # grabs density data        
    def density_row(row):
        if re.match('Density \(g', row[0]):
            densities.append(row[2])

    # grabs elemental data, splits into isotopes if need be        
    def elemental_row(row):
        if re.match('[A-Z][a-z]?-?(\d{1,3})?$', row[0]):
            element = nucname.id(row[0])
            weight_frac = row[3]
            if nucname.name(element) in elemental_mats:
                composition.update(elemental_mats[row[0]])
            else:
                composition[element] = float(weight_frac)
            nucids.add(element)
        else:
            pass

    # terminates collection of composition data, creates material        
    def ending_row(row):
        if re.match('Total$', row[0]):
            mat = Material(composition)
            mats.append(mat)
        else:
            pass

    # opens .csv, parses it
    with open(location, 'r') as f:
        reader = csv.reader(f)
        composition = {}
        name = ''
        for row in reader:
            starting_row(row)
            density_row(row)
            elemental_row(row)
            ending_row(row)

# Writes to file
def make_materials_compendium(nuc_data):
    """Adds materials compendium to nuc_data.h5."""
    # open nuc_data, make nuc_zz an array
    filters = tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)
    with tb.openFile(nuc_data, 'r+', filters=filters) as f:
        f.createGroup('/', 'material_library')
        f.createArray('/material_library', 'nucid', np.array(sorted(nucids)))

    # Writes elements for which we have compositional data to file
    for zz in elemental_mats:
        if 0 == len(elemental_mats[zz]):
            continue
        element = Material(elemental_mats[zz], mass=1.0, 
                           attrs={'name': nucname.name(zz)})
        element.write_hdf5(nuc_data, datapath="/material_library/materials", 
                           nucpath="/material_library/nucid", chunksize=70)

    # Writes materials from mats to file, and names them.
    for i in range(len(mats)):
        mats[i].mass = 1.0
        mats[i].density = float(densities[i])
        mats[i].attrs = {'name': names[i]}
        mats[i].write_hdf5(nuc_data, datapath="/material_library/materials", 
                           nucpath="/material_library/nucid", chunksize=70)
    
def make_materials_library(args):
    """Controller function for adding materials library."""
    nuc_data = args.nuc_data
    if os.path.exists(nuc_data):
        with tb.openFile(nuc_data, 'r') as f:
            if '/material_library' in f:
                print "skipping materials library data table creation; already exists."
                return

    # First make the elements
    print "Making the elements..."
    make_elements()

    # Then grab the materials compendium
    print "Grabbing materials compendium..."
    grab_materials_compendium(os.path.join(os.path.split(__file__)[0], 
                              'materials_compendium.csv'))

    # Make atomic weight table once we have the array
    print "Making materials library..."
    make_materials_compendium(nuc_data)
