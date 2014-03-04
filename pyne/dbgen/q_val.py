"""
Module allows the grabbing of q_values (energy per disintegration) for the calculation of decay heat. This currently consists of the nuclide, it's q_value, and the percent of q coming from gammas. This data is from 'ORIGEN-S DECAY DATA LIBRARY AND HALF-LIFE UNCERTAINTIES' (http://web.ornl.gov/~webworks/cppr/y2001/rpt/97914.pdf)
"""

from __future__ import print_function
import csv

import numpy as np
import tables as tb

from pyne import nucname

# Parses data from .csv
def grab_q_values(location = 'q_val_actinides.csv'):
    """Parses data from three q_val csv files."""

    # create list
    all_q_vals = []
    
    # grabs data row by row
    def read_row(row):
        if row[0] == 'Nuclide':
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
        all_q_vals.append(entry)
    
    # opens .csv files and parses them
    with open(location, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            read_row(row)
            
    return all_q_vals        

#def make_q_value_table(q_data):


def make_q_value():
    """Controller function for adding q-values"""
    q_values = args.q_values
    if os.path.exists(q_values):
        with tb.openFile(q_values, 'r') as f:
            if '/q_values' in f:
                print("skipping q_value table creation; already exists.")
                return

    # Grab the q_values
    print("Grabbing q_values...")
    grab_q_values(os.path.join(os.path.split(__file__)[0], 
                              'q_val_actinides.csv')) 

    # Make q_value table once we have the array
#    print("Making q_value table...")
#    make_q_value_table(q_values)
