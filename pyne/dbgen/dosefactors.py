"""
Module allows the grabbing of dose rate factors for the calculation of radiotoxicity. There are four dose rates provided: 
 (1) external from air (mrem/h per Ci/m^3)
     Table includes: nuclide, air dose rate factor, ratio to inhalation dose (All EPA values)
 (2) external from 15 cm of soil (mrem/h per Ci/m^2)
     Table includes: nuclide, GENII, EPA, DOE, GENII/EPA, DOE/EPA
 (3) ingestion (mrem/pCi)
     Table includes: nuclide, f1 (fraction of the activity ingested that enters body fluids.), GENII, EPA, DOE, GENII/EPA, DOE/EPA
 (4) inhalation (mrem/pCi)
     Table includes: nuclide, lung model*, GENII, EPA, DOE, GENII/EPA, DOE/EPA 

This data is from:

Appendix O 
[Exposure Scenarios and Unit Dose Factors for the Hanford Immobilized Low-Activity Tank Waste Performance Assessment, ref. HNF-SD-WM-TI-707 Rev. 1 December 1999]

of HNF-5636 
[DATA PACKAGES FOR THE HANFORDIMMOBILIZED LOW-ACTIVITY TANK WASTE PERFORMANCE ASSESSMENT: 2001 VERSION]

*Lung Model:
  "V" for tritium stands for vapor (50% larger absorption)
  "Organic" for C-14 means that the carbon is assumed to have an organic chemical form. 
  "D" material clears the lungs in days
  "W" material clears the lungs in weeks
  "Y" material clears the lungs in years
"""

from __future__ import print_function
import csv
import os

import numpy as np
import tables as tb

from pyne import nucname
from pyne.api import nuc_data
from pyne.dbgen.api import BASIC_FILTERS

# Parses data from .csv
def grab_dose_factors():
    """Parses data from dose factor csv files.
    """
    
    # Loops through four files
    df_superlist = []
    df_files = ['dosefactors_external_air.csv', 'dosefactors_external_soil.csv', 'dosefactors_ingest.csv', 'dosefactors_inhale.csv']
    for fname in df_files:
        dose_factors = []
        # Opens .csv file and parses it
        with open(os.path.join(os.path.dirname(__file__), fname), 'r') as f:
            reader = csv.reader(f)
            next(f)
            for row in reader:
                entry = read_row(row)
                dose_factors.append(entry)
        df_superlist.append(dose_factors)

    # Create four dose factor lists
    ext_air_df = df_superlist[0]
    ext_soil_df = df_superlist[1]
    ingest_df = df_superlist[2]
    inhale_df = df_superlist[3]
    
    return ext_air_df, ext_soil_df, ingest_df, inhale_df

# Grabs data row by row
def read_row(row):
    """Returns a list for each nuclide. Form varies based on type of dose rate factor:
    (1) External DF in Air: [int, float, float]
    (2) External DF in Soil: [int, float, float, float]
    (3) Ingestion DF: [int, float, float, float, float]
    (4) Inhalation DF: [int, string, float, float, float]
    
    Parameters
    ----------
    row : tuple
        One entry in a dose factor file.
    """

    # Create tuple
    entry = ()
    
    # Evaluate each component of the given row

    if row[0].endswith('+D'):
        row[0] = row[0][:-2]
    nuclide = nucname.id(row[0])

    # Case 1: DF from External Air
    if len(row) == 3:
        df_air = float(row[1])
        if len(row[2]) == 0:
            ratio = None
        else:
            ratio = float(row[2])
        entry = (nuclide, df_air, ratio)
    # Case 2: DF from External Soil
    elif len(row) == 6:
        genii = float(row[1])
        epa = float(row[2])
        doe = float(row[3])
        entry = (nuclide, genii, epa, doe)
    # Case 4: DF from Inhalation
    elif len(row) == 7 and row[1].isalpha():
        lungmodel = row[1]
        genii = float(row[2])
        epa = float(row[3])
        doe = float(row[4])
        entry = (nuclide, lungmodel, genii, epa, doe)
    # Case 4: DF from Inhalation
    else:
        f1 = float(row[1])
        genii = float(row[2])
        epa = float(row[3])
        doe = float(row[4])
        entry = (nuclide, f1, genii, epa, doe)
    
    return entry

# Write dose factor tables to file
def make_df_tables(ext_air_df, ext_soil_df, ingest_df, inhale_df, nuc_data, build_dir=""):
    """Adds four dose factor tables to the nuc_data.h5 library.

    Parameters
    ----------
    ext_air_df: list of tuples
        Array of external dose factors of air for tracked nuclides.
    ext_soil_df: list of tuples
        Array of external dose factors of soil for tracked nuclides.
    ingest_df: list of tuples
        Array of ingested dose factors for tracked nuclides.
    inhale_df: list of tuples
        Array of inhaled dose factors for tracked nuclides.
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place q_value files in.
    """
    
    # Define data types for all four cases
    ext_air_dtype = np.dtype([
        ('nuc', int),
        ('air_df', float),
        ('ratio', float),
        ])
    ext_soil_dtype = np.dtype([
        ('nuc', int),
        ('genii', float),
        ('epa', float),
        ('doe', float),
        ])
    ing_dtype = np.dtype([
        ('nuc', int),
        ('fluid_frac', float),
        ('genii', float),
        ('epa', float),
        ('doe', float),
        ])
    inh_dtype = np.dtype([
        ('nuc', int),
        ('lung_mod', 'S10'),
        ('genii', float),
        ('epa', float),
        ('doe', float),
        ])
    
    # Convert to numpy arrays
    ext_air_array = np.array(ext_air_df, dtype=ext_air_dtype)
    ext_soil_array = np.array(ext_soil_df, dtype=ext_soil_dtype)
    ingest_array = np.array(ingest_df, dtype=ing_dtype)
    inhale_array = np.array(inhale_df, dtype=inh_dtype)

    # Open the hdf5 file
    nuc_file = tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS)

    # Create a group for the tables
    df_group = nuc_file.createGroup("/neutron", "dose_factors", "Dose Rate Factors")

    # Make four new tables
    ext_air_table = nuc_file.createTable(df_group, 'external_air', ext_air_array, 'Nuclide, Air Dose Factor [mrem/h per Ci/m^3], Fraction of Air Dose Factor to Inhaled Dose Factor')
    ext_soil_table = nuc_file.createTable(df_group, 'external_soil', ext_soil_array, 'Nuclide, GENII [mrem/h per Ci/m^2], EPA, DOE')    
    ingest_table = nuc_file.createTable(df_group, 'ingestion', ingest_array, 'Nuclide, Frac of Activity in Body Fluids, GENII [mrem/pCi], EPA, DOE')    
    inhale_table = nuc_file.createTable(df_group, 'inhalation', inhale_array, 'Nuclide, Lung Model, GENII [mrem/pCi], EPA, DOE')

    # Ensure that data was written to table
    ext_air_table.flush()
    ext_soil_table.flush()
    ingest_table.flush()
    inhale_table.flush()

    # Close the hdf5 file
    nuc_file.close()

def make_dose_factors(args):
    """Controller function for adding dose factors"""
    nuc_data, build_dir = args.nuc_data, args.build_dir
    if os.path.exists(nuc_data):
        with tb.openFile(nuc_data, 'r') as f:
            if '/neutron/dose_factors' in f:
                print("skipping creation of dose factor tables; already exists.")
                return
    
    # Grab the dose factors from each file
    print('Grabbing dose factors...')
    ext_air_df, ext_soil_df, ingest_df, inhale_df = grab_dose_factors()
    
    # Make the four dose factor tables and writes them to file
    print("Making tables...")
    make_df_tables(ext_air_df, ext_soil_df, ingest_df, inhale_df, nuc_data, build_dir)
