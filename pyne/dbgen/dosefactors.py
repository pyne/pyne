"""
Module allows the grabbing of dose rate factors for the calculation of radiotoxicity. There are four dose rates provided:
 1. external from air (mrem/h per Ci/m^3)
     Table includes: nuclide, air dose rate factor, ratio to inhalation dose (All EPA values)
 2. external from 15 cm of soil (mrem/h per Ci/m^2)
     Table includes: nuclide, GENII, EPA, DOE, GENII/EPA, DOE/EPA
 3. ingestion (mrem/pCi)
     Table includes: nuclide, f1 (fraction of the activity ingested that enters body fluids.), GENII, EPA, DOE, GENII/EPA, DOE/EPA
 4. inhalation (mrem/pCi)
     Table includes: nuclide, lung model*, GENII, EPA, DOE, GENII/EPA, DOE/EPA 

This data is from: 
[Exposure Scenarios and Unit Dose Factors for the Hanford
Immobilized Low-Activity Tank Waste Performance Assessment, ref.
HNF-SD-WM-TI-707 Rev. 1 December 1999] Appendix O of HNF-5636 [DATA PACKAGES
FOR THE HANFORD IMMOBILIZED LOW-ACTIVITY TANK WASTE PERFORMANCE ASSESSMENT:
2001 VERSION]

Liability Disclaimer:
The PyNE Development Team shall not be liable for any loss or injury resulting
from decisions made with this data. 

*Lung Model:
  "V" for tritium stands for vapor (50% larger absorption)
  "O" for C-14 means that the carbon is assumed to have an Organic chemical form. 
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


def read_row(row):
    """Returns a list for each nuclide. Form varies based on type of dose rate factor:
    1. External DF in Air: [int, float, float]
    2. External DF in Soil: [int, float, float, float]
    3. Ingestion DF: [int, float, float, float, float]
    4. Inhalation DF: [int, string, float, float, float]

    Parameters
    ----------
    row : tuple
        One entry in a dose factor file.
    """

    # Create list
    entry = []

    # Evaluate each component of the given row
    if row[0].endswith("+D"):
        row[0] = row[0][:-2]
    nuclide = nucname.id(row[0])

    # Case 1: DF from External Air
    if len(row) == 3:
        dose_air = float(row[1])
        if len(row[2]) == 0:
            ratio = None
        else:
            ratio = float(row[2])
        entry = [nuclide, dose_air, ratio]
    # Case 2: DF from External Soil
    elif len(row) == 6:
        genii = float(row[1])
        epa = float(row[2])
        doe = float(row[3])
        entry = [nuclide, genii, epa, doe]
    # Case 4: DF from Inhalation
    elif len(row) == 7 and row[1].isalpha():
        lungmodel = row[1]
        genii = float(row[2])
        epa = float(row[3])
        doe = float(row[4])
        entry = [nuclide, genii, epa, doe, lungmodel]
    # Case 4: DF from Ingestion
    else:
        f1 = float(row[1])
        genii = float(row[2])
        epa = float(row[3])
        doe = float(row[4])
        entry = [nuclide, genii, epa, doe, f1]
    return entry


def grab_dose_factors():
    """Parses data from dose factor csv files."""

    # Populates Dose Factor list with initial set of nuclides: opens first .csv file and parses it
    dose_factors = []
    with open(
        os.path.join(os.path.dirname(__file__), "dosefactors_external_air.csv"), "r"
    ) as f:
        reader = csv.reader(f)
        next(f)
        next(f)
        for row in reader:
            entry = read_row(row)
            dose_factors.append(entry)

    # Loops through remaining three files to add other dose factors to each nuclide
    dose_files = [
        "dosefactors_external_soil.csv",
        "dosefactors_ingest.csv",
        "dosefactors_inhale.csv",
    ]
    for fname in dose_files:
        # Opens remaining .csv files and parses them
        with open(os.path.join(os.path.dirname(__file__), fname), "r") as f:
            reader = csv.reader(f)
            next(f)
            next(f)
            for row in reader:
                entry = read_row(row)
                # Adds info to nuclide's row
                for nuclide in dose_factors:
                    if entry[0] == nuclide[0]:
                        nuclide += entry[1 : len(entry)]

    # Create three dose factor lists with respect to source
    genii = []
    epa = []
    doe = []
    for nuclide in dose_factors:
        for i, val in enumerate(nuclide):
            if val is None:
                nuclide[i] = -1
        genii_row = (
            nuclide[0],
            -1,
            -1,
            nuclide[3],
            nuclide[6],
            nuclide[9],
            nuclide[10],
            nuclide[13],
        )
        genii.append(genii_row)
        epa_row = (
            nuclide[0],
            nuclide[1],
            nuclide[2],
            nuclide[4],
            nuclide[7],
            nuclide[9],
            nuclide[11],
            nuclide[13],
        )
        epa.append(epa_row)
        doe_row = (
            nuclide[0],
            -1,
            -1,
            nuclide[5],
            nuclide[8],
            nuclide[9],
            nuclide[12],
            nuclide[13],
        )
        doe.append(doe_row)
    return genii, epa, doe


def make_dose_tables(genii, epa, doe, nuc_data, build_dir=""):
    """Adds three dose factor tables to the nuc_data.h5 library.

    Parameters
    ----------
    genii: list of tuples
        Array of dose factors calculated by the code GENII.
    epa: list of tuples
        Array of dose factors calculated by the EPA.
    doe: list of tuples
        Array of dose factors calculated by the DOE.
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place q_value files in.
    """

    # Define data types for all three cases
    dose_dtype = np.dtype(
        [
            ("nuc", int),
            ("ext_air_dose", float),
            ("ratio", float),
            ("ext_soil_dose", float),
            ("ingest_dose", float),
            ("fluid_frac", float),
            ("inhale_dose", float),
            ("lung_mod", "S10"),
        ]
    )

    # Convert to numpy arrays
    genii_array = np.array(genii, dtype=dose_dtype)
    epa_array = np.array(epa, dtype=dose_dtype)
    doe_array = np.array(doe, dtype=dose_dtype)

    # Open the hdf5 file
    nuc_file = tb.open_file(nuc_data, "a", filters=BASIC_FILTERS)

    # Create a group for the tables
    dose_group = nuc_file.create_group("/", "dose_factors", "Dose Rate Factors")

    # Make three new tables
    genii_table = nuc_file.create_table(
        dose_group,
        "GENII",
        genii_array,
        "Nuclide, External Air Dose Factor [mrem/h per Ci/m^3], Fraction of Ext Air Dose to Inhalation Dose, External Soil Dose Factor [mrem/h per Ci/m^2], Ingestion Dose Factor [mrem/pCi], Fraction of Activity in Body Fluids, Inhalation Dose Factor [mrem/pCi], Lung Model Used",
    )
    epa_table = nuc_file.create_table(
        dose_group,
        "EPA",
        epa_array,
        "Nuclide, External Air Dose Factor [mrem/h per Ci/m^3], Fraction of Ext Air Dose to Inhalation Dose, External Soil Dose Factor [mrem/h per Ci/m^2], Ingestion Dose Factor [mrem/pCi], Fraction of Activity in Body Fluids, Inhalation Dose Factor [mrem/pCi], Lung Model Used",
    )
    doe_table = nuc_file.create_table(
        dose_group,
        "DOE",
        doe_array,
        "Nuclide, External Air Dose Factor [mrem/h per Ci/m^3], Fraction of Ext Air Dose to Inhalation Dose, External Soil Dose Factor [mrem/h per Ci/m^2], Ingestion Dose Factor [mrem/pCi], Fraction of Activity in Body Fluids, Inhalation Dose Factor [mrem/pCi], Lung Model Used",
    )

    # Ensure that data was written to table
    genii_table.flush()
    epa_table.flush()
    doe_table.flush()

    # Close the hdf5 file
    nuc_file.close()


def make_dose_factors(args):
    """Controller function for adding dose factors"""

    nuc_data, build_dir = args.nuc_data, args.build_dir
    if os.path.exists(nuc_data):
        with tb.open_file(nuc_data, "r") as f:
            if "/dose_factors" in f:
                print("skipping creation of dose factor tables; already exists.")
                return

    # Grab the dose factors from each file
    print("Grabbing dose factors...")
    genii, epa, doe = grab_dose_factors()

    # Make the 3 dose factor tables and writes them to file
    print("Making dose factor tables...")
    make_dose_tables(genii, epa, doe, nuc_data, build_dir)
