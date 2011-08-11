"""This module provides a way to grab and store raw data for atomic weights."""
import os
import re

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

    natural_nuclides = set()

    # Grab and parse elemental summary files.
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



def other_stuff():
    

    isotab = []

    for key in isolist:
        AW = 0.0
        AW_sig = 0.0
        Abund = 0.0			

        NucFetched = False

        while not NucFetched:
            try:
                print key
                kaeri = urllib2.urlopen( 'http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=%s'%(key) )
                NucFetched = True
            except:
                print "Failed to grab, retrying",

        for line in kaeri:
            if 0 < line.count("Atomic Mass:"):
                ls = line.split()
                AW = ls[2]
                AW_sig = ls[4]
            elif 0 < line.count("Atomic Percent Abundance:"):
                ls = line.split()
                abund_try = ls[-1]
                while Abund == 0.0:
                    try:
                        Abund = float(abund_try) / 100.0
                    except:
                        abund_try = abund_try[:-1]
        kaeri.close()

        if AW == 0.0:
            continue

        isotab.append([isoname.LLZZZM_2_aazzzm(key), AW, AW_sig, '%G'%Abund])


    isotab = sorted(isotab)

    libfile = open(file_out, 'w')
    for row in isotab:
        new_row = '{0:<6}  {1:<11}  {2:<9}  {3}\n'.format(isoname.aazzzm_2_LLZZZM(row[0]), row[1], row[2], row[3])
        libfile.write(new_row)
    libfile.close()

    return


"""Functions to make a nuclear data hdf5 file from the raw libraries."""
import os
import re
import math

import numpy as np
import tables as tb

import isoname

############################
### Next, Atomic Weights ###
############################

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
