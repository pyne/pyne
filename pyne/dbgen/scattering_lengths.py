"""This module provides a way to grab and store raw data for neutron scattering 
lengths.  This data comes from Neutron News, Vol. 3, No. 3, 1992, pp. 29-37 via 
a NIST webpage (http://www.ncnr.nist.gov/resources/n-lengths/list.html).  Please
contact Alan Munter, <alan.munter@nist.gov> for more information."""

import os
import re
import shutil
import urllib2

import numpy as np
import tables as tb

from pyne import nucname



def grab_scattering_lengths(build_dir="", file_out='scattering_lengths.html'):
    """Grabs the scattering cross-section lengths for neutrons from the NIST website
    or locally from this module."""
    build_filename = os.path.join(build_dir, file_out)
    local_filename = os.path.join(os.path.split(__file__)[0], file_out)

    if os.path.exists(local_filename):
        shutil.copy(local_filename, build_filename)
        return 

    nist = urllib2.urlopen("http://www.ncnr.nist.gov/resources/n-lengths/list.html")
    with open(build_filename, 'w') as f:
        f.write(nist.read())




def make_scattering_lengths(nuc_data, build_dir):
    # Grab the raw data
    print "Grabbing the scattering length data."
    grab_scattering_lengths(build_dir)
