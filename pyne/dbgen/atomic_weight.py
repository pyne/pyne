"""This module provides a way to grab and store raw data for atomic weights."""
import os
import re
import urllib2

import numpy as np
import tables as tb

from pyne import nucname
from pyne.dbgen.kaeri import grab_kaeri_nuclide

# Note that since ground state and meta-stable isotopes are of the same atomic weight, 
# the meta-stables have been discluded from the following data sets.

nat_iso_regex = re.compile('.*?/cgi-bin/nuclide[?]nuc=([A-Za-z]{1,2}\d{1,3}).*?[(].*?[)]')

def parse_for_natural_isotopes(htmlfile):
    """Parses an elemental html file, returning a set of naturally occuring isotopes."""
    nat_isos = set()
    with open(htmlfile, 'r') as f:
        for line in f:
            m = nat_iso_regex.search(line)
            if m is not None:
                nat_isos.add(nucname.zzaaam(m.group(1)))
    return nat_isos


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

    # Grab and parse elemental summary files.
    natural_nuclides = set()
    for element in nucname.name_zz.keys():
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
amdc_regex = re.compile('[ \d-]*? (\d{1,3})[ ]{1,4}(\d{1,3}) [A-Z][a-z]? .*? (\d{1,3}) ([ #.\d]{10,11}) ([ #.\d]{1,10})[ ]*?$')

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

        nuc = int(m.group(1) + m.group(2)) * 10
        mass = float(m.group(3)) + 1E-6 * float(m.group(4).strip().replace('#', ''))
        error = 1E-6 * float(m.group(5).strip().replace('#', ''))

        atomic_masses.append((nuc, mass, error))
        
    f.close()

    return atomic_masses    




atomic_weight_desc = {
    'iso_LL': tb.StringCol(itemsize=6, pos=0),
    'iso_zz': tb.IntCol(pos=1),
    'value':  tb.FloatCol(pos=2),
    'error':  tb.FloatCol(pos=3),
    'abund':  tb.FloatCol(pos=4),
    }

def _make_atomic_weight(h5_file='nuc_data.h5', data_file='atomic_weight.txt'):
    """Makes an atomic weight table and adds it to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the atomic weight text file to load data from.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Make a new the table
    Atable = kdb.createTable("/", "A", atomic_weight_desc, "Atomic Weight Data [amu]")
    nuc = Atable.row

    with open(data_file, 'r') as f:
        for line in f:
            ls = line.split()
            iso_LL = isoname.mixed_2_LLAAAM(ls[0])
            iso_zz = isoname.LLAAAM_2_zzaaam(iso_LL)

            nuc['iso_LL'] = iso_LL
            nuc['iso_zz'] = iso_zz
            nuc['value'] = float(ls[1])
            nuc['error'] = float(ls[2])
            nuc['abund'] = float(ls[3])

            # Insert nuclide to table
            nuc.append()

    # Ensure that data was written to table
    Atable.flush()

    # Close the hdf5 file
    kdb.close()




def make_atomic_weight(nuc_data, build_dir):
    # First grab the atomic abundance data
    print "Grabing the atomic abundance from KAERI"
    grab_kaeri_atomic_abund(build_dir)

    # Then grab mass data
    print "Grabing atomic mass data from AMDC"
    grab_atmoic_mass_adjustment(build_dir)


    parse_atmoic_mass_adjustment(build_dir)
