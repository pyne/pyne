"""This module provides a way to grab and store raw data for radioactive decay."""
import os
import re
import urllib
import urllib2
from zipfile import ZipFile

import numpy as np
import tables as tb

from pyne import nucname


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
    ensdf_zip = ['ensdf_1010_099.zip', 'ensdf_1010_199.zip', 'ensdf_1010_294.zip',]

    for f in ensdf_zip:
        fpath = os.path.join(build_dir, f)
        if f not in os.listdir(build_dir):
            print "  grabing {0} and placing it in {1}".format(f, fpath)
            urllib.urlretrieve(iaea_base_url + f, fpath)

        with ZipFile(fpath) as zf:
            for name in zf.namelist():
                print "    extracting {0} from {1}".format(name, f)
                zf.extract(name, build_dir)



#half_life_regex = re.compile('<li>Half life: [~]?(\d+[.]?\d*?)\s*(\w+)')

def parse_decay(build_dir=""):
    """Builds and returns a list of nuclide decay data.
    build_dir = os.path.join(build_dir, 'KAERI')

    # Grab and parse elemental summary files.
    nuclides = set()
    for element in nucname.name_zz.keys():
        htmlfile = element + '.html'
        nuclides = nuclides | parse_for_all_isotopes(os.path.join(build_dir, htmlfile))

    decay_data = []

    from_nuc_name = ""
    from_nuc_zz = 0
    to_nuc_name = ""
    to_nuc_zz = 0
    hl = 0.0
    dc = 0.0
    br = 1.0

    for nuc in nuclides:
        nuc_name = nucname.name(nuc)
        htmlfile = os.path.join(build_dir, nuc_name + '.html')

        from_nuc_name = nuc_name
        from_nuc_zz = nuc
        br = 1.0

        with open(htmlfile, 'r') as f:
            for line in f:
                m = half_life_regex.search(line)
                if m is not None:
                    val = float(m.group(1)) * 0.01
                    atomic_abund[nuc] = val
                    continue
    """
    # stub for now
    decay_data = [('H2', 10020, 'H1', 10010, 10.0, 0.1, 1.0)]
    return decay_data





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




def make_decay(nuc_data, build_dir):
    #with tb.openFile(nuc_data, 'r') as f:
    #    if hasattr(f.root, 'decay'):
    #        return 

    # grab the decay data
    print "Grabing the ENSDF decay data from IAEA"
    grab_ensdf_decay(build_dir)

    # Make atomic weight table once we have the array
    #print "Making atomic weight data table."
    #make_atomic_weight_table(nuc_data, build_dir)

