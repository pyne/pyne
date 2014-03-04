"""
Module allows the grabbing of q_values (energy per disintegration) for the calculation of decay heat. This currently consists of the nuclide, it's q_value, and the percent of q coming from gammas. This data is from 'ORIGEN-S DECAY DATA LIBRARY AND HALF-LIFE UNCERTAINTIES' (http://web.ornl.gov/~webworks/cppr/y2001/rpt/97914.pdf)
"""

from __future__ import print_function
import csv
import re

import numpy as np
import tables as tb

from pyne import nucname

# Parses data from .csv
def grab_q_value_data(location = 'q_val_actinides.csv'):
    """Parses data from three q_val csv files."""

    # create dictionary
    all_q_vals = {}
    
    # grabs data row by row
    def read_row(row):
        if re.match('Nuclide', row[0]):
            pass
        nuclide = nucname.id(join(row[0, 1]))
        q_val = row[2]
        gamma_frac = row[3]
        all_q_vals.append(nuclide, q_val, gamma_frac)
    
    # opens .csv files and parses them
    with open(location, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            read_row(row)
        

#def make_q_value_table(q_data):


def make_q_val():
    """Controller function for adding q-values"""
    q_data = args.q_data
    if os.path.exists(q_data):
        with tb.openFile(q_data, 'r') as f:
            if '/q_val_data' in f:
                print("skipping q_value data table creation; already exists.")
                return

    # Grab the q_value data
    print("Grabbing q_value data...")
    grab_q_value_data(os.path.join(os.path.split(__file__)[0], 
                              'q_val_actinides.csv')) 

    # Make q_value table once we have the array
#    print("Making q_value table...")
#    make_q_value_table(q_data)
