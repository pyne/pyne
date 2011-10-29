"""This module provides a way to grab and store raw data for atomic weights."""
import os
import re
import urllib2

import numpy as np
import tables as tb

from pyne import nucname
from pyne.dbgen.kaeri import grab_kaeri_nuclide, parse_for_natural_isotopes

# Note that since ground state and meta-stable isotopes are of the same atomic weight, 
# the meta-stables have been discluded from the following data sets.


def grab_kaeri_atomic_abund(build_dir=""):
    """Grabs the KAERI files needed for the atomic abundance calculation, 
    if not already present.

    Parameters
    ----------
    build_dir : str
        Major directory to place html files in. 'KAERI/' will be appended.
    """
    # Add kaeri to build_dir
    build_dir = os.path.join(build_dir, 'KAERI')
    try:
        os.makedirs(build_dir)
    except OSError:
        pass
    already_grabbed = set(os.listdir(build_dir))

    # Grab and parse elemental summary files.
    natural_nuclides = set()
    for element in nucname.name_zz.keys():
        htmlfile = element + '.html'
        if htmlfile not in already_grabbed:
            grab_kaeri_nuclide(element, build_dir)

        natural_nuclides = natural_nuclides | parse_for_natural_isotopes(os.path.join(build_dir, htmlfile))

    # Grab natural nuclide files
    for nuc in natural_nuclides:
        nuc = nucname.name(nuc)
        htmlfile = nuc + '.html'
        if htmlfile not in already_grabbed:
            grab_kaeri_nuclide(nuc, build_dir)



atomic_abund_regex = re.compile('<li>Atomic Percent Abundance: (\d+[.]?\d*?)%')

def parse_atomic_abund(build_dir=""):
    """Builds and returns a dictionary from nuclides to atomic abundence fractions."""
    build_dir = os.path.join(build_dir, 'KAERI')

    # Grab and parse elemental summary files.
    natural_nuclides = set()
    for element in nucname.name_zz.keys():
        htmlfile = element + '.html'
        natural_nuclides = natural_nuclides | parse_for_natural_isotopes(os.path.join(build_dir, htmlfile))

    atomic_abund = {}    

    for nuc in natural_nuclides:
        nuc_name = nucname.name(nuc)
        htmlfile = os.path.join(build_dir, nuc_name + '.html')

        with open(htmlfile, 'r') as f:
            for line in f:
                m = atomic_abund_regex.search(line)
                if m is not None:
                    val = float(m.group(1)) * 0.01
                    atomic_abund[nuc] = val
                    break

    return atomic_abund



def grab_atmoic_mass_adjustment(build_dir=""):
    """Grabs the current atomic mass adjustment from the Atomic
    Mass Data Center.  These are courtesy of Georges Audi and 
    Wang Meng via a private communication, April 2011."""
    mass_file = 'mass.mas114'
    bd_files = os.listdir(build_dir)
    if mass_file in bd_files:
        return 

    mass = urllib2.urlopen('http://amdc.in2p3.fr/masstables/Ame2011int/mass.mas114')
    with open(os.path.join(build_dir, mass_file), 'w') as f:
        f.write(mass.read())


# Note, this regex specifically leaves our free neutrons
#amdc_regex = re.compile('[ \d-]*? (\d{1,3})[ ]{1,4}(\d{1,3}) [A-Z][a-z]? .*? (\d{1,3}) ([ #.\d]{10,11}) ([ #.\d]{1,10})[ ]*?$')
amdc_regex = re.compile('[ \d-]*? (\d{1,3})[ ]{1,4}(\d{1,3}) [A-Z][a-z]? .*? (\d{1,3}) ([ #.\d]{5,12}) ([ #.\d]+)[ ]*?$')

def parse_atmoic_mass_adjustment(build_dir=""):
    """Parses the atomic mass adjustment data into a list of tuples of 
    the nuclide, atomic mass, and error."""
    mass_file = 'mass.mas114'
    f = open(os.path.join(build_dir, mass_file), 'r')

    atomic_masses = []

    for line in f:
        m = amdc_regex.search(line)
        if m is None:
            continue

        nuc = (10000 * int(m.group(1))) + (10 * int(m.group(2)))
        mass = float(m.group(3)) + 1E-6 * float(m.group(4).strip().replace('#', ''))
        error = 1E-6 * float(m.group(5).strip().replace('#', ''))

        atomic_masses.append((nuc, mass, error))
        
    f.close()

    return atomic_masses    




atomic_weight_desc = {
    'nuc_name': tb.StringCol(itemsize=6, pos=0),
    'nuc_zz':   tb.IntCol(pos=1),
    'mass':     tb.FloatCol(pos=2),
    'error':    tb.FloatCol(pos=3),
    'abund':    tb.FloatCol(pos=4),
    }

atomic_weight_dtype = np.dtype([
    ('nuc_name', 'S6'),
    ('nuc_zz',   int),
    ('mass',     float),
    ('error',    float), 
    ('abund',    float), 
    ])

def make_atomic_weight_table(nuc_data, build_dir=""):
    """Makes an atomic weight table in the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place html files in.
    """
    # Grab raw data
    atomic_abund  = parse_atomic_abund(build_dir)
    atomic_masses = parse_atmoic_mass_adjustment(build_dir)

    A = {}

    # Add normal isotopes to A
    for nuc_zz, mass, error in atomic_masses:
        try: 
            nuc_name = nucname.name(nuc_zz)
        except RuntimeError:
            continue

        if nuc_zz in atomic_abund:
            A[nuc_zz] = nuc_name, nuc_zz, mass, error, atomic_abund[nuc_zz]
        else:
            A[nuc_zz] = nuc_name, nuc_zz, mass, error, 0.0

    # Add naturally occuring elements
    for element in nucname.name_zz:
        nuc_zz = nucname.zzaaam(element)
        A[nuc_zz] = element, nuc_zz, 0.0, 0.0, 0.0
        
    for nuc, abund in atomic_abund.items():
        zz = nuc / 10000
        element_zz = zz * 10000
        element = nucname.zz_name[zz]

        nuc_name, nuc_zz, nuc_mass, _error, _abund = A[nuc]
        elem_name, elem_zz, elem_mass, _error, _abund = A[element_zz]

        new_elem_mass = elem_mass + (nuc_mass * abund)
        A[element_zz] = element, element_zz, new_elem_mass, 0.0, 0.0


    A = sorted(A.values(), key=lambda x: x[1])
    #A = np.array(A, dtype=atomic_weight_dtype)

    # Open the HDF5 File
    kdb = tb.openFile(nuc_data, 'a')

    # Make a new the table
    Atable = kdb.createTable("/", "atomic_weight", atomic_weight_desc, 
                             "Atomic Weight Data [amu]", expectedrows=len(A))
    Atable.append(A)

    # Ensure that data was written to table
    Atable.flush()

    # Close the hdf5 file
    kdb.close()




def make_atomic_weight(nuc_data, build_dir):
    """Controller function for adding atomic_weights."""
    if os.path.exists(nuc_data):
        with tb.openFile(nuc_data, 'r') as f:
            if hasattr(f.root, 'atomic_weight'):
                return 

    # First grab the atomic abundance data
    print "Grabing the atomic abundance from KAERI"
    grab_kaeri_atomic_abund(build_dir)

    # Then grab mass data
    print "Grabing atomic mass data from AMDC"
    grab_atmoic_mass_adjustment(build_dir)

    # Make atomic weight table once we have the array
    print "Making atomic weight data table."
    make_atomic_weight_table(nuc_data, build_dir)

