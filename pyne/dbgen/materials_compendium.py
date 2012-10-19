import csv
import re
from pyne.material import Material
import tables as tb
import numpy as np
from pyne import nucname
import os
import pyne.data as pd

nuc_zz = set()
mats = []
names = []
elems = set()


def starting_row(row):
	if re.match('\d{1,3}\.  ', row[0]):
		#print row
		name = row[1]
		names.append(name)
		composition.clear()
		#print name
	else:
		pass
		
def elemental_row(row):
	if re.match('[A-Z][a-z]?-?(\d{1,3})?$', row[0]):
		#print row
		element = nucname.zzaaam(row[0])
		weight_frac = row[3]
		composition[element] = float(weight_frac)
		nuc_zz.add(element)
	else:
		pass


		
def ending_row(row):
	if re.match('Total$', row[0]):
		#print name, composition, '\n', len(composition)
		mat = Material(composition)
		mats.append(mat)
		#print mat, "\n"
	else:
		pass

with open('materials_compendium.csv', 'r') as f:
    reader = csv.reader(f)
    composition = {}
    name = ''
    for row in reader:
        starting_row(row)
        elemental_row(row)
        ending_row(row)

with tb.openFile('test.h5', 'w', filters=tb.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=False)) as f:
	f.createArray('/', 'nuc_zz', np.array(list(nuc_zz)))

for nuc in nuc_zz:
	elem = nuc/10000
	abund = pd.natural_abund(nuc)
	elems.add(elem)

for i in range(len(mats)):
	mats[i].attrs = {'name': names[i]}
	mats[i].write_hdf5('test.h5', chunksize = 62) #12:258.3s
	#print mats[i], "\n"
