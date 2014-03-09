"""
Module allows the grabbing of q_values (energy per disintegration) for the calculation of decay heat. This currently consists of the nuclide, it's q_value, and the percent of q coming from gammas. This data is from 'ORIGEN-S DECAY DATA LIBRARY AND HALF-LIFE UNCERTAINTIES' (http://web.ornl.gov/~webworks/cppr/y2001/rpt/97914.pdf)
"""

from __future__ import print_function
import csv
import os

import numpy as np
import tables as tb

from pyne import nucname
from pyne.api import nuc_data
#from pyne.api import BASIC_FILTERS

# Parses data from .csv
def grab_q_values(fname):
    """Parses data from three q_val csv files.
    
    Parameters
    ----------
    fname : str
        Path to q_value file.
    """
    
    # Create list
    all_q_values = []
        
    # Grabs data row by row
    def read_row(row):
        if row[0] == 'Nuclide' or len(row[0].strip()) == 0:
            return
        nuclide = nucname.id(''.join(row[0:2]).replace(' ', ''))
        if len(row[2]) == 0:
            q_val = 0.0
        else:
            q_val = float(row[2])
        if len(row[3]) == 0:
            gamma_frac = 0.0
        else:   
            gamma_frac = float(row[3])
        entry = [nuclide, q_val, gamma_frac]
        all_q_values.append(entry)
    
    # Opens .csv files and parses them
    with open(fname, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            read_row(row)

    return all_q_values         

# Sorts and filters list of q_values
def format_q_values(all_q_values):
    """Filters the q_value data for multiple entries then sorts the nuclides.
    
    Parameters
    ----------
    all_q_values : list of lists
        Array of q_values for all nuclides.
    """

    distinct_all_q_values = []
    
    # Ensures only one entry per nuclide
    for nuclide in all_q_values:
        if not nuclide in distinct_all_q_values:
            distinct_all_q_values.append(nuclide)
    
    # Sort in order of nuclide
    distinct_all_q_values.sort(key=lambda nucid: nucid[0])
 
    return distinct_all_q_values

# Write q_value table to file
def make_q_value_table(all_q_values, nuc_data, build_dir=""):
    """Adds q_value table to the nuc_data.h5 library.

    Parameters
    ----------
    all_q_values: list of lists
        Array of q_values for all nuclides.
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place q_value files in.
    """

    # Sort and filter the q_values
    distinct_all_q_values = format_q_values(all_q_values)
    
    # Iterate over all nuclides in the q_value table
    for row in distinct_all_q_values:

        row = (nuclide['Nuclide'],
               nuclide['Q_value [MeV per disintegration]'],
               nuclide['Fraction of Q that comes from gammas'])

        distinct_all_q_values.append(row)

    # Converts to numpy array
    q_value_array = np.array(distinct_all_q_values)

    # Open the hdf5 file and create group
    nuc_file = tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS)
    q_val_group = nuc_file.createGroup('/', 'q_values', 'Q_values for nuclides')

    # Make a new table
    q_value_table = nuc_file.createTable(q_val_group, np.empty(0), 
                    'Nuclide, Q_value [MeV per disintegration],' 
                    'Fraction of Q that comes from gammas')
    q_value_table.append(q_value_array)

    # Ensure that data was written to table
    q_value_table.flush()

    # Close the hdf5 file
    q_values.close()

def make_q_value(args):
    """Controller function for adding q-values"""
#    q_values = args.q_values
#    if os.path.exists(q_values):
#        with tb.openFile(q_values, 'r') as f:
#            if '/q_values' in f:
#                print("skipping q_value table creation; already exists.")
#                return

    nuc_data, build_dir = args.nuc_data, args.build_dir
    
    # Grab the q_values
    print('Grabbing q_values...')
    q_value_files = ['q_val_actinides.csv', 'q_val_fissionproducts.csv', 
                     'q_val_light.csv']
    all_q_values = []
    for fname in q_value_files:
        all_q_values += grab_q_values(fname)
    
    # Make the q_value table and write to file
    print("Making q_value table...")
    make_q_value_table(all_q_values, nuc_data, build_dir))     
