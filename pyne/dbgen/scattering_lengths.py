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




def _init_scattering_length(kdb):
    """Initializes the scattering length part of the database.

    Keyword Args:
        * kdb (tables.File): a nuclear data hdf5 file.
    """

    # Create neutron group
    if not hasattr(kdb.root, 'neutron'):
        neutron_group = kdb.createGroup('/', 'neutron', 'Neutron Cross Sections')


nist_iso_pattern = "(\d*)([A-Za-z]+)"

def nist_2_zzaaam(nist_iso):
    """Converts a NIST style isotope to a zzaaam style.

    Args:
        * nist_iso (str): A nist isotope.

    Returns:
        * iso_zz (int): a zzaaam isotope.
    """
    m = re.match(nist_iso_pattern, nist_iso)

    if m.group(1) == "":
        elem = m.group(2).upper()
        iso_zz = isoname.LLzz[elem] * 10000
    else:
        iso_zz = isoname.mixed_2_zzaaam(m.group(2) + m.group(1))

    return iso_zz


def nist_2_LLAAAM(nist_iso):
    """Converts a NIST style isotope to a LLAAAM style.

    Args:
        * nist_iso (str): A nist isotope.

    Returns:
        * iso_LL (int): a LLAAAM isotope.
    """
    m = re.match(nist_iso_pattern, nist_iso)

    iso_LL = isoname.mixed_2_LLAAAM(m.group(2) + m.group(1))

    return iso_LL


def nist_num(nist_data):
    """Converts a NIST style data point to a point.

    Args:
        * nist_data (str): A nist data point.

    Returns:
        * d (int): a data point.
    """
    nd = nist_data
    while ('(' in nd) or (')' in nd):
        nd_pre = nd.partition('(')
        nd_post = nd.partition(')')
        nd = nd_pre[0] + nd_post[2]

    if (nd == "---") or (nd == ""):
        nd = "0.0"

    nd = nd.replace("<i>i</i>", 'j')
    d = eval(nd)
    return d

scat_len_desc = {
    'iso_LL': tb.StringCol(6, pos=0),
    'iso_zz': tb.Int32Col(pos=1),

    'b_coherent': tb.ComplexCol(16, pos=2),
    'b_incoherent': tb.ComplexCol(16, pos=3),

    'xs_coherent': tb.Float64Col(pos=4),
    'xs_incoherent': tb.Float64Col(pos=5),
    'xs': tb.Float64Col(pos=6),
    }

scat_len_data = "[ aeE()<>i/.+\d-]+?"
scat_len_space = "[ \t]+"
scat_len_pattern = "<td>{space}(?P<iso>[A-Za-z\d]+){space}<td>{space}(?P<conc>{data}){space}<td>{space}(?P<b_coherent>{data}){space}<td>{space}(?P<b_incoherent>{data}){space}<td>{space}(?P<xs_coherent>{data}){space}<td>{space}(?P<xs_incoherent>{data}){space}<td>{space}(?P<xs>{data}){space}<td>{space}(?P<xs_a>{data}){space}<tr>".format(data=scat_len_data, space=scat_len_space)

def make_scattering_lengths(h5_file='nuc_data.h5', data_file='scattering_lengths.html'):
    """Adds the neutron fission product yields to the hdf5 library.

    Keyword Args:
        * h5_file (str): path to hdf5 file.
        * data_file (str): path to the scattering_length.html data file.
    """
    # Open the HDF5 File
    kdb = tb.openFile(h5_file, 'a')

    # Ensure that the appropriate file structure is present
    _init_scattering_length(kdb)

    # Read in cinder data file
    with open(data_file, 'r') as f:
        raw_data = f.read()

    # Init the neutron fission product info table
    sl_table = kdb.createTable('/neutron/', 'scattering_lengths', scat_len_desc, 
                                'Neutron Scattering Lengths, b [cm], sigma [barns]')
    slrow = sl_table.row

    # Iterate over all isotopes in the table
    for m in re.finditer(scat_len_pattern, raw_data):
        md = m.groupdict()

        # Adds the row to the table
        slrow['iso_LL'] = nist_2_LLAAAM(md['iso'])
        slrow['iso_zz'] = nist_2_zzaaam(md['iso'])

        slrow['b_coherent'] = nist_num(md['b_coherent']) * (10**-13)
        slrow['b_incoherent'] = nist_num(md['b_incoherent']) * (10**-13)

        slrow['xs_coherent'] = nist_num(md['xs_coherent'])
        slrow['xs_incoherent'] = nist_num(md['xs_incoherent'])
        slrow['xs'] = nist_num(md['xs'])

        slrow.append()

    # Write the table
    sl_table.flush()

    # Close the hdf5 file
    kdb.close()



def make_scattering_lengths(nuc_data, build_dir):
    # Grab the raw data
    print "Grabbing the scattering length data."
    grab_scattering_lengths(build_dir)


