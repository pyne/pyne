"""Module handles parsing EAF formatted cross section files and adding
the data to PyNE's HDF5 storage.
"""

import os
import re
import shutil
from glob import glob

import numpy as np
import tables as tb

from pyne import nucname
from pyne.utils import to_barns, failure
from pyne.dbgen.api import BASIC_FILTERS


# numpy array row storage information for EAF data
eaf_dtype = np.dtype([
    ('nuc_zz',        int          ),
    ('rxnum',         'S7'         ),
    ('rxstr',         'S4'         ),
    ('daughter',      'S5'         ),
    ('xsec',          float, (175,))
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
    eaf_file = "/filespace/groups/cnerg/opt/FENDL2.0-A/fendlg-2.0_175"

    with open(eaf_file, 'r') as f:
        raw_data = f.read()

    eaf_data = list()
    eaf_pattern = eaf_info_pattern + eaf_bin_pattern 

    # Iterate over all iso/rx combinations in file
    for m in re.finditer(eaf_pattern, raw_data, re.DOTALL):
        md = m.groupdict()

        xsec_list = [float(x) for x in md['xsec'].split()]
        xsec_list += (175-len(xsec_list))*[0.0]

        # Store information in new row of array.
        eafrow = (
                  nucname.zzaaam(md['iso']),
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

    # Add group structure by calling placeholder function to get boundaries
    db.createArray('/neutron/eaf_xs', 'E_g', _get_eaf_groups(), \
            'Neutron energy group bounds [MeV]')
    
    # Close the HDF5 file
    db.close()


def _get_eaf_groups():
    """Function for retrieving EAF group structure.

    Note, energies are entered high to low in the HDF5 storage.

    This is used as a placeholder for a future PyNE implementation of a database
    storing group structures.  This hypothetical database would return the
    group structure like this method does.

    Returns
    ----------
    eaf_E_g_array : 1D list
        List of energy group boundaries from high to low.
    
    """

    eaf_E_g_array = [0.0,
            1.0000E-07, 4.1399E-07, 5.3158E-07, 6.8256E-07, 8.7643E-07,
            1.1254E-06, 1.4450E-06, 1.8554E-06, 2.3824E-06, 3.0590E-06,
            3.9279E-06, 5.0435E-06, 6.4760E-06, 8.3153E-06, 1.0677E-05,
            1.3710E-05, 1.7604E-05, 2.2603E-05, 2.9023E-05, 3.7267E-05,
            4.7851E-05, 6.1442E-05, 7.8893E-05, 1.0130E-04, 1.3007E-04,
            1.6702E-04, 2.1445E-04, 2.7536E-04, 3.5358E-04, 4.5400E-04,
            5.8295E-04, 7.4852E-04, 9.6112E-04, 1.2341E-03, 1.5846E-03,
            2.0347E-03, 2.2487E-03, 2.4852E-03, 2.6126E-03, 2.7465E-03,
            3.0354E-03, 3.3546E-03, 3.7074E-03, 4.3074E-03, 5.5308E-03,
            7.1017E-03, 9.1188E-03, 1.0595E-02, 1.1709E-02, 1.5034E-02,
            1.9305E-02, 2.1875E-02, 2.3579E-02, 2.4176E-02, 2.4788E-02,
            2.6058E-02, 2.7000E-02, 2.8501E-02, 3.1828E-02, 3.4307E-02,
            4.0868E-02, 4.6309E-02, 5.2475E-02, 5.6562E-02, 6.7380E-02,
            7.2025E-02, 7.9499E-02, 8.2503E-02, 8.6517E-02, 9.8037E-02,
            1.1109E-01, 1.1679E-01, 1.2277E-01, 1.2907E-01, 1.3569E-01,
            1.4264E-01, 1.4996E-01, 1.5764E-01, 1.6573E-01, 1.7422E-01,
            1.8316E-01, 1.9255E-01, 2.0242E-01, 2.1280E-01, 2.2371E-01,
            2.3518E-01, 2.4724E-01, 2.7324E-01, 2.8725E-01, 2.9452E-01,
            2.9721E-01, 2.9849E-01, 3.0197E-01, 3.3373E-01, 3.6883E-01,
            3.8774E-01, 4.0762E-01, 4.5049E-01, 4.9787E-01, 5.2340E-01,
            5.5023E-01, 5.7844E-01, 6.0810E-01, 6.3928E-01, 6.7206E-01,
            7.0651E-01, 7.4274E-01, 7.8082E-01, 8.2085E-01, 8.6294E-01,
            9.0718E-01, 9.6167E-01, 1.0026E+00, 1.1080E+00, 1.1648E+00,
            1.2246E+00, 1.2874E+00, 1.3534E+00, 1.4227E+00, 1.4957E+00,
            1.5724E+00, 1.6530E+00, 1.7377E+00, 1.8268E+00, 1.9205E+00,
            2.0190E+00, 2.1225E+00, 2.2313E+00, 2.3069E+00, 2.3457E+00,
            2.3653E+00, 2.3851E+00, 2.4660E+00, 2.5924E+00, 2.7253E+00,
            2.8651E+00, 3.0119E+00, 3.1664E+00, 3.3287E+00, 3.6788E+00,
            4.0657E+00, 4.4933E+00, 4.7237E+00, 4.9659E+00, 5.2205E+00,
            5.4881E+00, 5.7695E+00, 6.0653E+00, 6.3763E+00, 6.5924E+00,
            6.7032E+00, 7.0469E+00, 7.4082E+00, 7.7880E+00, 8.1873E+00,
            8.6071E+00, 9.0484E+00, 9.5123E+00, 1.0000E+01, 1.0513E+01,
            1.1052E+01, 1.1618E+01, 1.2214E+01, 1.2523E+01, 1.2840E+01,
            1.3499E+01, 1.3840E+01, 1.4191E+01, 1.4550E+01, 1.4918E+01,
            1.5683E+01, 1.6487E+01, 1.6905E+01, 1.7333E+01, 1.9640E+01]
    
    return eaf_E_g_array


def make_eaf(args):
    """Controller function for adding cross section data from EAF format file.
    """
    nuc_data, build_dir, datapath = args.nuc_data, args.build_dir, args.datapath

    # Check if the table already exists
    with tb.openFile(nuc_data, 'a', filters=BASIC_FILTERS) as f:
        if hasattr(f.root, 'neutron') and hasattr(f.root.neutron, 'eaf_xs'):
            return

    #
    print "Grabbing the EAF activation data."
    parse_eaf_xsec(build_dir)

    print "Making EAF activation data table."
    make_eaf_table(nuc_data, build_dir)

