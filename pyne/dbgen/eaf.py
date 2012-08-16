""" """
import os
import re
import shutil
from glob import glob

import numpy as np
import tables as tb

from pyne import nucname
from pyne.utils import to_barns, failure
from pyne.dbgen.api import BASIC_FILTERS

<<<<<<< HEAD
=======
def grab_eaf_data(build_dir=""):
    """Grabs the EAF file from the local filesystem if not already present."""
    build_filename = os.path.join(build_dir, 'eaf')
    if os.path.exists(build_filename):
        return True

    # FIXME make local_filename more general
    local_filename = "/filespace/groups/cnerg/opt/FENDL2.0-A/fendlg-2.0_175"
    if not os.path.exists(local_filename):
        print failure("EAF file not found - skipping.")
        return False

    print "Grabbing the EAF activation data from " + local_filename
    shutil.copy(local_filename, build_filename)
    return True


#def _init_eaf(db):
#    """Initializes a multigroup cross-section part of the database.
#
#    Parameters
#    ----------
#    db : tables.File 
#        A nuclear data hdf5 file.
#    """
#
#    # Create neutron group
#    if not hasattr(db.root, 'neutron'):
#        neutron_group = db.createGroup('/', 'neutron', 'Neutron Interaction Data')
#
#    # Create xs group
#    if not hasattr(db.root.neutron, 'eaf_xs'):
#        nxs_mg_group = db.createGroup("/neutron", "eaf_xs", "EAF 175-Group Neutron Activation Cross Section Data")
#
#    # Create fission_yield groups
#    if not hasattr(db.root.neutron, 'cinder_fission_products'):
#        nxs_mg_group = db.createGroup("/neutron", "cinder_fission_products", "CINDER Neutron Fission Product Yield Data")


>>>>>>> pub-origin/staging

# numpy array row storage information for EAF data
eaf_dtype = np.dtype([
<<<<<<< HEAD
    ('nuc_zz',        int          ),
    ('rxnum',         'S7'         ),
    ('rxstr',         'S4'         ),
    ('daughter',      'S5'         ),
    ('xsec',          float, (175,))
=======
    ('nuc', int),
    ('rxnum', 'S7'),
    ('rxstr', 'S4'),
    ('daughter', 'S5'),
    ('xsec', float, (175,))
>>>>>>> pub-origin/staging
    ])

# Regular expression for parsing an individual set of EAF data
# Includes some groupnames that are currently unused.
eaf_info_pattern = \
    "(?P<iso>\d{5,7})\s*(?P<rxnum>\d{2,4})\s*(?P<ngrps>\d{1,3})" \
    + "\s*(?P<parent>[a-zA-Z]{1,2}\s{0,3}\d{1,3}[M ][12 ])" \
    + "(?P<rxstr>\(N,[\w\s]{3}\))(?P<daugh>[a-zA-Z.]{1,2}\s{0,3}\d{1,3})(.*?)"
eaf_bin_pattern = "(?P<xsec>(\d\.\d{5}E[-+]\d{2}\s*){1,175})"


def parse_eaf_xsec(build_dir):
    """Create numpy array by parsing EAF data

    Parameters
    ----------
    build_dir : str
        Directory where EAF data is stored (TODO).
    
    Returns
    ---------
    eaf_array : numpy array
        Numpy array with a row for each isotope+reaction combination
         found in the EAF data.

    """
    
    #TODO: change eaf_file to something more universal
<<<<<<< HEAD
    eaf_file = "/filespace/groups/cnerg/opt/FENDL2.0-A/fendlg-2.0_175"
=======
    eaf_file = os.path.join(build_dir, 'eaf')
>>>>>>> pub-origin/staging

    with open(eaf_file, 'r') as f:
        raw_data = f.read()

    eaf_data = list()
    eaf_pattern = eaf_info_pattern + eaf_bin_pattern 

    # Iterate over all iso/rx combinations in file
    for m in re.finditer(eaf_pattern, raw_data, re.DOTALL):
        md = m.groupdict()

        xsec_list = [float(x) for x in md['xsec'].split()]
        xsec_list += (175-len(xsec_list))*[0.0]

<<<<<<< HEAD
        # Store information in new row of array.
        eafrow = (
                  nucname.zzaaam(md['iso']),
=======
        eafrow = (nucname.zzaaam(md['iso']),
>>>>>>> pub-origin/staging
                  md['rxnum'],
                  md['rxstr'],
                  md['daugh'],
                  xsec_list
                  )
        
        eaf_data.append(eafrow)

    eaf_array = np.array(eaf_data, dtype=eaf_dtype)

    print "Read in {0} sets of EAF data.".format(len(eaf_array))

    return eaf_array


def make_eaf_table(nuc_data, build_dir=""):
    """Function for adding EAF group and table to HDF5 storage.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory where EAF data is located.
    
    """
    eaf_array = parse_eaf_xsec(build_dir)

    # Open the HDF5 file
    db = tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS)

    # Ensure that the appropriate file structure is present
    if not hasattr(db.root, 'neutron'):
        neutron_group = db.createGroup('/', 'neutron', \
                'Neutron Interaction Data')

    # Create xs group
    if not hasattr(db.root.neutron, 'eaf_xs'):
        eaf_group = db.createGroup("/neutron", "eaf_xs", \
            "EAF 175-Group Neutron Activation Cross Section Data")

    eaf_table = db.createTable("/neutron/eaf_xs", "eaf_xs", \
            np.empty(0, dtype=eaf_dtype), \
            "EAF Activation Cross Section Data [barns]", \
            expectedrows=len(eaf_array))

    # Put eaf_array in the table
    eaf_table.append(eaf_array)

    # Write the table
    eaf_table.flush()

    # Close the HDF5 file
    db.close()


def make_eaf(args):
<<<<<<< HEAD
    """Controller function for adding cross section data from EAF format file.
    """
    nuc_data, build_dir, datapath = args.nuc_data, args.build_dir, args.datapath
=======
    """Controller function for adding cinder data."""
    nuc_data, build_dir = args.nuc_data, args.build_dir
>>>>>>> pub-origin/staging

    # Check if the table already exists
    with tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS) as f:
        if hasattr(f.root, 'neutron') and hasattr(f.root.neutron, 'eaf_xs'):
            return

<<<<<<< HEAD
    #
    print "Grabbing the EAF activation data."
    parse_eaf_xsec(build_dir)
=======
    # First grab the EAF activation data.
    grabbed = grab_eaf_data(build_dir)
    if not grabbed:
        return
>>>>>>> pub-origin/staging

    print "Making EAF activation data table."
    make_eaf_table(nuc_data, build_dir)
