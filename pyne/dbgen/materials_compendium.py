import csv	
import re
from pyne.material import Material
import tables as tb
import numpy as np
from pyne import nucname
import os
from pyne.data import natural_abund, natural_abund_map
habund = natural_abund('H')

nuc_zz = set()
mats = []
names = []
densities = []

# Make a dictionary that represents elements as materials composed of
# their isotopes
elemental_mats = {}
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

# opens file, parses it
with open('materials_compendium.csv', 'r') as f:
    reader = csv.reader(f)
    composition = {}
    name = ''
    for row in reader:
        starting_row(row)
        density_row(row)
        elemental_row(row)
        ending_row(row)

# make nuc_zz an array
with tb.openFile('test.h5', 'w', filters=tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)) as f:
	f.createArray('/', 'nuc_zz', np.array(list(nuc_zz)))

# Writes elements for which we have compositional data to file
for zz in elemental_mats:
	if 0 == len(elemental_mats[zz]):
		continue
	element = Material(elemental_mats[zz], mass = 1.0, attrs = {'name':nucname.name(zz)})
	element.write_hdf5('test.h5', chunksize = 70)

# Writes materials from mats to file, and names them.
for i in range(len(mats)):
	mats[i].mass = 1.0
	mats[i].density = float(densities[i])
	mats[i].attrs = {'name': names[i]}
	mats[i].write_hdf5('test.h5', chunksize = 70)
