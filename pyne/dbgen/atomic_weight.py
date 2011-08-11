"""This module provides a way to grab and store raw data for atomic weights."""

import isoname
import urllib2

from urllib import urlopen
from urllib import urlencode


# Note that since ground state and meta-stable isotopes are of the same atomic weight, 
# the meta-stables have been discluded from the following data sets.

def grab_kaeri_atomic_weights(file_out='atomic_weight.txt'):
    """Makes the atomic weight library.
    Library rows have the the following form:

    iso	AW	AW_sig	Abund

    where:
        iso	= Isotope in LLZZZM format
        AW	= Atomic Weight [amu]
        AW_sig	= Atomic Weight Uncertainty [amu]
        Abund	= Natural fractional atomic abundance [unitless]

    Not to be used under normal circumstances.
    More like an embedded script, in case the librrary file is lost and unrecoverable.

    FIXME: This could use a rewrite such that it doesn't have to grab them all at once.
    """

    isolist = []
    
    for key in isoname.LLaadic.keys():
        NucFetched = False

        while not NucFetched:
            try:
                print key 
                kaeri = urllib2.urlopen( 'http://atom.kaeri.re.kr/cgi-bin/nuclide?nuc=%s'%(key) )
                NucFetched = True
            except:
                print "Failed to grab, retrying",

        for line in kaeri:
            if 0 < line.count("/cgi-bin/nuclide?nuc="):
                nuc = line.partition("/cgi-bin/nuclide?nuc=")[2].partition("\"")[0].upper()
                if not (nuc == key):
                    isolist.append(nuc)
        kaeri.close()

    print "\n~~~~~~~~\n"

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

def make_atomic_weight(h5_file='nuc_data.h5', data_file='atomic_weight.txt'):
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


