"""This module provides a way to grab and store raw data for radioactive decay."""
import os
import re
import glob
import urllib
import urllib2
from zipfile import ZipFile

import numpy as np
import tables as tb

from pyne import nucname
from pyne import ensdf
from pyne.dbgen.api import BASIC_FILTERS

def grab_ensdf_decay(build_dir=""):
    """Grabs the ENSDF decay data files
    if not already present.

    Parameters
    ----------
    build_dir : str
        Major directory to place html files in. 'KAERI/' will be appended.
    """
    # Add kaeri to build_dir
    build_dir = os.path.join(build_dir, 'ENSDF')
    try:
        os.makedirs(build_dir)
    except OSError:
        pass

    # Grab ENSDF files and unzip them.
    iaea_base_url = 'http://www-nds.iaea.org/ensdf_base_files/2010-November/'
    s3_base_url = 'http://s3.amazonaws.com/pyne/'
    ensdf_zip = ['ensdf_1010_099.zip', 'ensdf_1010_199.zip', 'ensdf_1010_294.zip',]

    for f in ensdf_zip:
        fpath = os.path.join(build_dir, f)
        if f not in os.listdir(build_dir):
            print "  grabbing {0} and placing it in {1}".format(f, fpath)
            urllib.urlretrieve(iaea_base_url + f, fpath)

            if os.path.getsize(fpath) < 1048576: 
                print "  could not get {0} from IAEA; trying S3 mirror".format(f)
                os.remove(fpath)
                urllib.urlretrieve(s3_base_url + f, fpath)

        # not using ZipFile context manager (with statement for Python 2.6)
        try:
            zf = ZipFile(fpath)
            for name in zf.namelist():
                if not os.path.exists(os.path.join(build_dir, name)):
                    print "    extracting {0} from {1}".format(name, f)
                    zf.extract(name, build_dir)
        finally:
            zf.close()
        



atomic_decay_dtype = np.dtype([
    ('from_nuc',   int),
    ('level',         float),
    ('to_nuc',     int), 
    ('half_life',     float),
    ('decay_const',   float),
    ('branch_ratio',  float),
    ])
    
    
def parse_decay(build_dir=""):
    """Builds and returns a list of nuclide decay data.
    """
    build_dir = os.path.join(build_dir, 'ENSDF')

    decay_data = []
    files = sorted([f for f in glob.glob(os.path.join(build_dir, 'ensdf.*'))])
    for f in files:
        print "    parsing decay data from {0}".format(f)
        decay_data += ensdf.half_life(f)

    ln2 = np.log(2.0)
    decay_data = [(fn, lvl, tn, hl, ln2/hl, br) for fn, lvl, tn, hl, br in decay_data]
    decay_data = set(decay_data)
    decay_data = sorted(decay_data, key=lambda x: (x[1], x[4]))

    decay_array = np.array(decay_data, dtype=atomic_decay_dtype)
    #da, mask = np.unique(decay_array, return_index=True)
    #mask.sort()
    #decay_array = decay_array[mask]
    return decay_array
    
    
def make_atomic_decay_table(nuc_data, build_dir=""):
    """Makes a decay table in the nuc_data library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory to place ensdf files in.
    """
    # Grab raw data
    atomic_decay  = parse_decay(build_dir)

    # Open the HDF5 File
    decay_db = tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS)

    # Make a new the table
    decaytable = decay_db.createTable("/", "atomic_decay", 
                    np.empty(0, dtype=atomic_decay_dtype), 
                    "Atomic Decay Data level [MeV], half_life [s], decay_const "
                    "[1/s], branch_ratio [frac]", expectedrows=len(atomic_decay))
    decaytable.append(atomic_decay)

    # Ensure that data was written to table
    decaytable.flush()

    # Close the hdf5 file
    decay_db.close()




def make_decay(args):
    """Controller function for adding decay data."""
    nuc_data, build_dir = args.nuc_data, args.build_dir

    with tb.openFile(nuc_data, 'r') as f:
        if hasattr(f.root, 'atomic_decay'):
            return 

    # grab the decay data
    print "Grabbing the ENSDF decay data from IAEA"
    grab_ensdf_decay(build_dir)

    # Make atomic weight table once we have the array
    print "Making decay data table."
    make_atomic_decay_table(nuc_data, build_dir)

