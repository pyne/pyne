#!/usr/bin/env python

import urllib2

def get_atomic_masses():
    """Creates a dictionary of atomic masses based off of data from G Audi,
    A. H. Wapstra, and C. Thibault, "The 2003 Atomic Mass Evaluation
    (II). Tables, graphs, and references." Nuclear Physics A729, 337 (2003).
    
    Returns
    -------
    atomic_mass : dict
        A dictionary where each key is an integer equal to the ZAID identifier
        for the nuclide and the value is the atomic mass in units of amu
        (1.660538921e−27 kg). For example, to get the atomic mass of U-235 which
        has Z=92, you would enter atomic_mass[92235].
    """

    # Open table of masses from NNDC
    url = 'http://www.nndc.bnl.gov/masses/mass.mas03'
    mass_file = urllib2.urlopen(url)

    # Skip first 39 lines to get to data
    n_skip = 39
    [mass_file.readline() for i in range(n_skip)]

    # Create empty dictionary to store atomic mass
    atomic_mass = {}

    # Read data from file
    for line in mass_file:
        # Read atomic number and mass number
        Z = int(line[9:14])
        A = int(line[14:19])
        zaid = Z*1000 + A

        # Read atomic mass of this nuclide. Note that the mass is split into two
        # fields, one which has the integer part and the other which has the
        # fractional part.
        mass_int_part = float(line[96:99])
        mass_frac_part = float(line[99:111].replace('#','.'))

        # Calculate mass and add to dictionary
        mass = mass_int_part + 1e-6*mass_frac_part
        atomic_mass[zaid] = mass

    # Check to make sure mass of C-12 is 12 amu.
    assert atomic_mass[6012] == 12.0

    return atomic_mass
