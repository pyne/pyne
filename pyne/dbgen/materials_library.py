import csv	
import re
from pyne.material import Material
import tables as tb
import numpy as np
from pyne import nucname
import os
from pyne.data import natural_abund, natural_abund_map

elemental_mats = {}

# Make a dictionary that represents elements as dicts of their isotopes
def make_elements():
	habund = natural_abund('H')
	for name, zz in nucname.name_zz.items():
		elemental_mats[name] = {}
	for nuc, abund in natural_abund_map.items():
		if 0 == nuc%10000 or abund == 0.0:
			continue
		zz = nuc/10000
		if zz not in nucname.zz_name:
			continue
		name = nucname.zz_name[zz]
		elemental_mats[name][nuc] = abund
    
# Parses data from .csv
def grab_materials_compendium(location = 'materials_compendium.csv'):
	nuc_zz = set()
	mats = []
	names = []
	densities = []	
	# grabs name from starting row, starts a new dictionary for composition
	def starting_row(row):
		if re.match('\d{1,3}\.  ', row[0]):
			#print row
			name = row[1]
			names.append(name)
			composition.clear()
			#print name
		else:
			pass
	# grabs density data		
	def density_row(row):
		if re.match('Density \(g', row[0]):
			densities.append(row[2])
	# grabs elemental data, splits into isotopes if need be		
	def elemental_row(row):
		if re.match('[A-Z][a-z]?-?(\d{1,3})?$', row[0]):
			element = nucname.zzaaam(row[0])
			weight_frac = row[3]
			if nucname.name(element) in elemental_mats:
				composition.update(elemental_mats[row[0].upper()])
			else:
				composition[element] = float(weight_frac)
			nuc_zz.add(element)
		else:
			pass
	# terminates collection of composition data, creates material		
	def ending_row(row):
		if re.match('Total$', row[0]):
			#print name, composition, '\n', len(composition)
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
	# open nuc_data, make nuc_zz an array
	#with tb.openFile(nuc_data, 'r+', filters=tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)) as f:
	with tb.openFile(nuc_data, 'a') as f:
		f.createGroup('/', 'materials_library')
		f.createArray('/materials_library', 'nuc_zz', np.array(list(nuc_zz)))

	# Writes elements for which we have compositional data to file
	for zz in elemental_mats:
		if 0 == len(elemental_mats[zz]):
			continue
		element = Material(elemental_mats[zz], mass = 1.0, attrs = {'name':nucname.name(zz)})
		element.write_hdf5(nuc_data, datapath="/materials_library/materials", nucpath="/materials_library/nuc_zz", chunksize=70)

	# Writes materials from mats to file, and names them.
	for i in range(len(mats)):
		mats[i].mass = 1.0
		mats[i].density = float(densities[i])
		mats[i].attrs = {'name': names[i]}
		mats[i].write_hdf5(nuc_data, datapath="/material_library/materials", nucpath="/material_library/nuc_zz", chunksize=70)
	
def make_materials_library(args):
    """Controller function for adding materials library."""
    nuc_data = args.nuc_data

    if os.path.exists(nuc_data):
        with tb.openFile(nuc_data, 'r') as f:
            if '/materials_library' in f:
                return
        f.close()

    # First make the elements
    print "Making the elements"
    elemental_mats = {}
    make_elements()

    # Then grab the materials compendium
    print "Grabbing materials compendium"
    grab_materials_compendium(os.path.join(os.path.split(__file__)[0], 'materials_compendium.csv'))

    # Make atomic weight table once we have the array
    print "Making materials compendium"
    make_materials_compendium(nuc_data)
