"""This module provides a way to locate, parse, and store CINDER cross sections."""
import os
import re
import shutil
from glob import glob

import numpy as np
import tables as tb

from pyne import nucname
from pyne.utils import to_barns


def grab_cinder_dat(build_dir=""):
    """Grabs the cinder.dat file from the DATAPATH directory if not already present."""
    build_filename = os.path.join(build_dir, 'cinder.dat')
    if os.path.exists(build_filename):
        return 

    if 'DATAPATH' in os.environ:
        datapath = os.environ['DATAPATH']
        print "Grabing cinder.dat from " + datapath
    else:
        raise OSError("DATAPATH not defined in environment; cinder.dat not found.")

    local_filename = os.path.join(datapath, "[Cc][Ii][Nn][Dd][Ee][Rr].[Dd][Aa][Tt]")
    local_filename = glob(local_filename)
    if 0 < len(local_filename):
        shutil.copy(local_filename[0], build_filename)
    else:
        raise OSError("cinder.dat file not found in DATAPATH dir.")


# These read in cinder.dat
cinder_float = "[\d.+-Ee]+"

def _init_cinder(db):
    """Initializes a multigroup cross-section part of the database.

    Parameters
    ----------
    db : tables.File 
        A nuclear data hdf5 file.
    """

    # Create neutron and photon groups
    if not hasattr(db.root, 'neutron'):
        neutron_group = db.createGroup('/', 'neutron', 'Neutron Interaction Data')

    if not hasattr(db.root, 'photon'):
        photon_group = db.createGroup('/', 'photon', 'Photon Interaction Data')

    # Create xs group
    if not hasattr(db.root.neutron, 'cinder_xs'):
        nxs_mg_group = db.createGroup("/neutron", "cinder_xs", "CINDER Multi-Group Neutron Cross Section Data")

    # Create source groups
    if not hasattr(db.root.photon, 'cinder_source'):
        gxs_mg_group = db.createGroup("/photon", "cinder_source", "CINDER Multi-Group Photon Source Data")

    # Create fission_yield groups
    if not hasattr(db.root.neutron, 'cinder_fission_products'):
        nxs_mg_group = db.createGroup("/neutron", "cinder_fission_products", "CINDER Neutron Fission Product Yield Data")

    if not hasattr(db.root.photon, 'cinder_fission_products'):
        nxs_mg_group = db.createGroup("/photon", "cinder_fission_products", "CINDER Photofission Product Yield Data")




def get_group_sizes(raw_data):
    """Gets the number of nuclides and groups in this data file.

    Parameters
    ----------
    raw_data : str 
        Input cinder.dat data file as a string.

    Returns
    -------
    nuclides : int 
        The number of nuclides in the dataset
    G_n : int 
        The number of neutron energy groups in the dataset
    G_p : int 
        The number of proton energy groups in the dataset
    G_g : int 
        The number of photon energy groups in the dataset
    """
    # Search for the group pattern
    G_pattern = "(\d+) nuclides,\s+(\d+) neutron groups,\s+(\d+) proton groups,\s+(\d+) photon groups"
    m = re.search(G_pattern, raw_data)
    g = m.groups()

    # Convert to ints
    nuclides = int(g[0])
    G_n = int(g[1])
    G_p = int(g[2])
    G_g = int(g[3])

    return nuclides, G_n, G_p, G_g



def make_mg_group_structure(nuc_data, build_dir=""):
    """Add the energy group bounds arrays to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = get_group_sizes(raw_data)

    # Find & write neutron group structure
    n_E_g_pattern = "Neutron group .*, MeV" + ("\s+("+cinder_float+")")*(G_n + 1)
    m = re.search(n_E_g_pattern, raw_data)
    g = m.groups()
    n_E_g = np.array(g[::-1], dtype=float)
    db.createArray('/neutron/cinder_xs', 'E_g', n_E_g, 'Neutron energy group bounds [MeV]')

    # Find & write photon group structure
    g_E_g_pattern = "Gamma structure, MeV" + ("\s+("+cinder_float+")")*(G_g + 1)
    m = re.search(g_E_g_pattern, raw_data)
    g = m.groups()
    g_E_g = np.array(g[::-1], dtype=float)
    db.createArray('/photon/cinder_source', 'E_g', g_E_g, 'Photon energy group bounds [MeV]')

    # Close the hdf5 file
    db.close()


# Helpful patterns
from_nuc_pattern = "\n#[\s\d]{4}:\s+(\d+).*?\n(_______________________| [\w/-]+ Fission Yield Data)"
to_nuc_base = "#[\s\d]{4}:\s+(\d+) produced by the following C-X  \((.{4})\) REF:.*?\n"

absorption_dtype_tuple = [
    ('from_nuc_name', 'S6'),
    ('from_nuc_zz', int),
    ('to_nuc_name', 'S6'),
    ('to_nuc_zz', int),
    ('reaction_type', 'S4'),
    # Extra 'xs' entry should be appended with value ('xs', float, G_n)    
    ]

def make_mg_absorption(nuc_data, build_dir=""):
    """Adds the absorption reaction rate cross sections to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = get_group_sizes(raw_data)

    # Init the neutron absorption table
    absorption_dtype = np.dtype(absorption_dtype_tuple + [('xs', float, G_n)])
    absorption_table = db.createTable('/neutron/cinder_xs/', 'absorption', absorption_dtype, 
                                       'Neutron absorption reaction cross sections [barns]')
    abrow = absorption_table.row

    # Init to_nuc_pattern
    to_nuc_pattern = to_nuc_base + ("\s+("+cinder_float+")")*G_n

    # Iterate through all from nuctopes.
    for m_from in re.finditer(from_nuc_pattern, raw_data, re.DOTALL):
        from_nuc_zz = nucname.zzaaam(m_from.group(1))

        # Check matestable state
        if 1 < from_nuc_zz%10:
            # Metastable state too high!
            continue
        from_nuc_name = nucname.name(from_nuc_zz)

        # Grab the string for this from_nuc in order to get all of the to_nucs
        from_nuc_part = m_from.group(0)

        # Iterate over all to_nucs
        for m_to in re.finditer(to_nuc_pattern, from_nuc_part):
            to_nuc_zz = nucname.zzaaam(m_to.group(1))

            # Check matestable state
            if 1 < to_nuc_zz%10:
                # Metastable state too high!
                continue
            to_nuc_name = nucname.name(to_nuc_zz)

            # Munge reaction type
            rx_type = m_to.group(2)
            rx_type = rx_type.strip()

            # Setup XS array
            xs = np.array(m_to.groups()[-1:1:-1], dtype=float)
            assert xs.shape == (G_n, )

            # Write this row to the absorption table
            abrow['from_nuc_name'] = from_nuc_name
            abrow['from_nuc_zz'] = from_nuc_zz

            abrow['to_nuc_name'] = to_nuc_name
            abrow['to_nuc_zz'] = to_nuc_zz

            abrow['reaction_type'] = rx_type
            abrow['xs'] = xs

            abrow.append()

        # Flush this from nuc
        absorption_table.flush()

    # Close the hdf5 file
    db.close()


fission_base = "\n([\s\d]+),([\s\d]+),([\s\d]+)=\(n,f\) yield sets\. If > 0,[\s\d]{4}-gp fisn CX follows\..*?\n"

fission_dtype_tuple = [
    ('nuc_name', 'S6'),
    ('nuc_zz', int),
    ('thermal_yield', np.int8),
    ('fast_yield', np.int8),
    ('high_energy_yield', np.int8),
    # Extra 'xs' entry should be appended with value ('xs', float, G_n)    
    ]

def make_mg_fission(nuc_data, build_dir=""):
    """Adds the fission reaction rate cross sections to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = get_group_sizes(raw_data)

    # Init the neutron absorption table
    fission_dtype = np.dtype(fission_dtype_tuple + [('xs', float, G_n)])
    fission_table = db.createTable('/neutron/cinder_xs/', 'fission', fission_dtype, 
                                    'Neutron fission reaction cross sections [barns]')
    frow = fission_table.row

    # Init to_nuc_pattern
    fission_pattern = fission_base + ("\s+("+cinder_float+")")*G_n

    # Iterate through all from nuctopes.
    for m_from in re.finditer(from_nuc_pattern, raw_data, re.DOTALL):
        from_nuc_zz = nucname.zzaaam(m_from.group(1))

        # Check matestable state
        if 1 < from_nuc_zz%10:
            # Metastable state too high!
            continue
        from_nuc_name = nucname.name(from_nuc_zz)

        # Grab the string for this from_nuc in order to get all of the to_nucs
        from_nuc_part = m_from.group(0)

        # Grab the fission part
        m_fission = re.search(fission_pattern, from_nuc_part)
        if m_fission is None:
            continue

        # Grab yield indexes
        yield_t = int(m_fission.group(1))
        yield_f = int(m_fission.group(2))
        yield_h = int(m_fission.group(3))

        # Grab XS array
        xs = np.array(m_fission.groups()[-1:2:-1], dtype=float)
        assert xs.shape == (G_n, )

        # Write fission table row
        frow['nuc_name'] = from_nuc_name
        frow['nuc_zz'] = from_nuc_zz

        frow['thermal_yield']     = yield_t
        frow['fast_yield']        = yield_f
        frow['high_energy_yield'] = yield_h

        frow['xs'] = xs

        # Write out this row
        frow.append()
        fission_table.flush()

    # Close the hdf5 file
    db.close()


gamma_decay_base = "gamma spectra from .{3} multiplied by\s+(" + cinder_float + ")\s+to agree w/ Eg=\s*(" +\
                   cinder_float + ")\n"

gamma_decay_dtype_tuple = [
    ('nuc_name', 'S6'),
    ('nuc_zz', int),
    ('energy', float),
    ('scaling_factor', float),
    # Extra 'spectrum' entry should be appended with value ('spectrum', float, G_g)
    ]

def make_mg_gamma_decay(nuc_data, build_dir=""):
    """Adds the gamma decay spectrum information to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    nuclides, G_n, G_p, G_g = get_group_sizes(raw_data)

    # Init the gamma absorption table
    gamma_decay_dtype = np.dtype(gamma_decay_dtype_tuple + [('spectrum', float, G_g)])
    gamma_decay_table = db.createTable('/photon/cinder_source/', 'decay_spectra', gamma_decay_dtype, 
                                        'Gamma decay spectrum [MeV]')
    gdrow = gamma_decay_table.row

    # Init to_nuc_pattern
    gamma_decay_pattern = gamma_decay_base + ("\s+("+cinder_float+")")*G_g

    # Iterate through all from nuctopes.
    for m_from in re.finditer(from_nuc_pattern, raw_data, re.DOTALL):
        from_nuc_zz = nucname.zzaaam(m_from.group(1))

        # Check matestable state
        if 1 < from_nuc_zz%10:
            # Metastable state too high!
            continue
        from_nuc_name = nucname.name(from_nuc_zz)

        # Grab the string for this from_nuc in order to get all of the to_nucs
        from_nuc_part = m_from.group(0)

        # Grab the fission part
        m_gd = re.search(gamma_decay_pattern, from_nuc_part)
        if m_gd is None:
            continue

        # Grab base data
        scale = float(m_gd.group(1))
        energy = float(m_gd.group(2))

        # Grab spectrum
        spectrum = np.array(m_gd.groups()[-1:1:-1], dtype=float)
        assert spectrum.shape == (G_g, )

        # Prepare the row
        gdrow['nuc_name'] = from_nuc_name
        gdrow['nuc_zz'] = from_nuc_zz

        gdrow['energy'] = energy
        gdrow['scaling_factor'] = scale

        gdrow['spectrum'] = spectrum

        # Write out the row
        gdrow.append()
        gamma_decay_table.flush()

    # Close the hdf5 file
    db.close()



def get_fp_sizes(raw_data):
    """Gets the number of fission product yield data sets in this file.

    Parameters
    ----------
    data : str
        Input cinder.dat data file as a string.

    Returns
    -------
    N_n : int 
        The number of neutron fission product yield datasets in the file.
    N_g int 
        The number of photon fission product yield datasets in the file.
    """
    # Search for the neutron pattern
    N_n_pattern = "Fission Yield Data.*?\n\s*(\d+)\s+yield sets"
    m_n = re.search(N_n_pattern, raw_data)
    N_n = int(m_n.group(1))

    # Search for the photon pattern
    N_g_pattern = "Photofission Yield Data.*?\n\s*(\d+)\s+yield sets"
    m_g = re.search(N_g_pattern, raw_data)
    N_g = int(m_g.group(1))

    return N_n, N_g


fp_info_dtype = np.dtype([
    ('index', np.int16),
    ('nuc_name', 'S6'),
    ('nuc_zz', int),
    ('type', 'S11'),
    ('mass', float),
    ])

fp_type_flag = {
    't': 'thermal', 
    'f': 'fast',
    'h': 'high_energy', 
    's': 'spontaneous', 
    }

iit_pattern = "(\d{1,3})\s+\d{2,3}-\s?([A-Z]{1,2}-[ Mm\d]{3})([tfhs])"
mass_pattern = "\d{1,3}\.\d{1,4}"

nfp_info_pattern = "Fission Yield Data.*?fission products"

def parse_neutron_fp_info(raw_data):
    """Grabs the neutron fission product info.

    Parameters
    ----------
    raw_data : str 
        string of the cinder.dat data file.

    Returns
    -------
    info_table : array 
        Structured array with the form "(index, nuc_name, nuc_zz, type, mass)". 
    """
    # Get group sizes
    N_n, N_g = get_fp_sizes(raw_data)

    # Grab the part of the file that is a neutron fission product yield info
    m_info = re.search(nfp_info_pattern, raw_data, re.DOTALL)
    nfp_info_raw = m_info.group(0)

    # Grab the index, nuctope, and type
    iits = re.findall(iit_pattern, nfp_info_raw)

    # Grab the masses 
    masses = re.findall(mass_pattern, nfp_info_raw)

    # Make sure data is the right size
    assert N_n == len(iits) 
    assert N_n == len(masses)

    # Make info table rows 
    info_table = []
    for m in range(N_n):
        iit = iits[m]
        index = int(iit[0])

        nuc_zz = nucname.zzaaam(iit[1])
        # Correct for metastable flag
        if 0 != nuc_zz%10:
            nuc_zz = nuc_zz + 2000

        nuc_name = nucname.name(nuc_zz)
        type = fp_type_flag[iit[2]]
        mass = float(masses[m])

        info_row = (index, nuc_name, nuc_zz, type, mass)
        info_table.append(info_row)

    info_table = np.array(info_table, dtype=fp_info_dtype)

    return info_table


def make_neutron_fp_info(nuc_data, build_dir=""):
    """Adds the neutron fission product yield info to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # get the info table
    info_table = parse_neutron_fp_info(raw_data)

    # Init the neutron fission product info table
    nfp_table = db.createTable('/neutron/cinder_fission_products/', 'info', fp_info_dtype, 
                                'CINDER Neutron Fission Product Yield Information')

    # Append Rows
    nfp_table.append(info_table)

    # Close the hdf5 file
    db.close()


gfp_info_pattern = "Photofission Yield Data.*?fission products"

def grab_photon_fp_info(raw_data):
    """Grabs the photon fission product info.

    Parameters
    ----------
    raw_data : str
        string of the cinder.dat data file.

    Returns
    -------
    info_table : array 
        Structured array with the form "(index, nuc_name, nuc_zz, type, mass)". 
    """
    # Get group sizes
    N_n, N_g = get_fp_sizes(raw_data)

    # Grab the part of the file that is a neutron fission product yield info
    m_info = re.search(gfp_info_pattern, raw_data, re.DOTALL)
    gfp_info_raw = m_info.group(0)

    # Grab the index, nuctope, and type
    iits = re.findall(iit_pattern, gfp_info_raw)

    # Grab the masses 
    masses = re.findall(mass_pattern, gfp_info_raw)

    # Make sure data is the right size
    assert N_g == len(iits) 
    assert N_g == len(masses)

    # Make info table rows 
    info_table = []
    for m in range(N_g):
        iit = iits[m]
        index = int(iit[0])

        nuc_zz = nucname.zzaaam(iit[1])
        # Correct for metastable flag
        if 0 != nuc_zz%10:
            nuc_zz = nuc_zz + 2000

        nuc_name = nucname.name(nuc_zz)
        type = fp_type_flag[iit[2]]
        mass = float(masses[m])

        info_row = (index, nuc_name, nuc_zz, type, mass)
        info_table.append(info_row)

    info_table = np.array(info_table, dtype=fp_info_dtype)

    return info_table

def make_photon_fp_info(nuc_data, build_dir=""):
    """Adds the photofission product yield info to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # Grab photon info table
    info_table = grab_photon_fp_info(raw_data)

    # Init the neutron fission product info table
    gfp_table = db.createTable('/photon/cinder_fission_products/', 'info', fp_info_dtype, 
                                'CINDER Photofission Product Yield Information')
    gfp_table.append(info_table)

    # Close the hdf5 file
    db.close()


fp_yields_dtype = np.dtype([
    ('index', np.int16),
    ('from_nuc_name', 'S6'),
    ('from_nuc_zz', int),
    ('to_nuc_name', 'S6'),
    ('to_nuc_zz', int),
    ('mass_frac', float),
    ])


fp_to_nuc_insert = "\s{1,3}" + cinder_float
fp_to_nuc_base = "  ([ \d]{4}) ([ \d]{7})\n("

nfp_yields_pattern = "Fission Yield Data.*Photofission Yield Data"

def make_neutron_fp_yields(nuc_data, build_dir=""):
    """Adds the neutron fission product yields to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    N_n, N_g = get_fp_sizes(raw_data)

    # get the info table
    info_table = parse_neutron_fp_info(raw_data)

    # Grab the part of the file that is a neutron fission product yields
    m_yields = re.search(nfp_yields_pattern, raw_data, re.DOTALL)
    nfp_yields_raw = m_yields.group(0)

    # Init the neutron fission product info table
    nfp_table = db.createTable('/neutron/cinder_fission_products/', 'yields', fp_yields_dtype, 
                                'CINDER Neutron Fission Product Yields')
    nfprow = nfp_table.row

    # Iterate over all to-nucs
    fp_to_nuc_pattern = fp_to_nuc_base + N_n*fp_to_nuc_insert + ")"
    for m_to in re.finditer(fp_to_nuc_pattern, nfp_yields_raw):
        to_nuc_zz = nucname.zzaaam(m_to.group(2).strip())

        # Check matestable state
        if 1 < to_nuc_zz%10:
            # Metastable state too high!
            continue
        to_nuc_name = nucname.name(to_nuc_zz)

        # Get the array of yield data
        yields = np.array(m_to.group(3).split(), dtype=float)
        assert len(yields) == N_n

        # Prep rows to the table
        for n in range(N_n):
            info = info_table[n]

            nfprow['index'] = info[0]
            nfprow['from_nuc_name'] = info[1]
            nfprow['from_nuc_zz'] = info[2]
            nfprow['to_nuc_name'] = to_nuc_name
            nfprow['to_nuc_zz'] = to_nuc_zz
            nfprow['mass_frac'] = yields[n]

            nfprow.append()

        # Write the table
        nfp_table.flush()

    # Close the hdf5 file
    db.close()


gfp_yields_pattern = "Photofission Yield Data.*"

def make_photon_fp_yields(nuc_data, build_dir):
    """Adds the photofission product yields to the hdf5 library.

    Parameters
    ----------
    nuc_data : str
        Path to nuclide data file.
    build_dir : str
        Directory with cinder.dat file.
    """
    # Open the HDF5 File
    db = tb.openFile(nuc_data, 'a')

    # Ensure that the appropriate file structure is present
    _init_cinder(db)

    # Read in cinder data file
    cinder_dat = os.path.join(build_dir, 'cinder.dat')
    with open(cinder_dat, 'r') as f:
        raw_data = f.read()

    # Get group sizes
    N_n, N_g = get_fp_sizes(raw_data)

    # get the info table
    info_table = grab_photon_fp_info(raw_data)

    # Grab the part of the file that is a neutron fission product yields
    m_yields = re.search(gfp_yields_pattern, raw_data, re.DOTALL)
    gfp_yields_raw = m_yields.group(0)

    # Init the neutron fission product info table
    gfp_table = db.createTable('/photon/cinder_fission_products/', 'yields', fp_yields_dtype, 
                                'CINDER Photofission Product Yields')
    gfprow = gfp_table.row

    # Iterate over all to-nucs
    fp_to_nuc_pattern = fp_to_nuc_base + N_g*fp_to_nuc_insert + ")"
    for m_to in re.finditer(fp_to_nuc_pattern, gfp_yields_raw):
        to_nuc_zz = nucname.zzaaam(m_to.group(2).strip())

        # Check matestable state
        if 1 < to_nuc_zz%10:
            # Metastable state too high!
            continue
        to_nuc_name = nucname.name(to_nuc_zz)

        # Get the array of yield data
        yields = np.array(m_to.group(3).split(), dtype=float)
        assert len(yields) == N_g

        # Prep rows to the table
        for n in range(N_g):
            info = info_table[n]

            gfprow['index'] = info[0]
            gfprow['from_nuc_name'] = info[1]
            gfprow['from_nuc_zz'] = info[2]
            gfprow['to_nuc_name'] = to_nuc_name
            gfprow['to_nuc_zz'] = to_nuc_zz
            gfprow['mass_frac'] = yields[n]

            gfprow.append()

        # Write the table
        gfp_table.flush()

    # Close the hdf5 file
    db.close()


def make_cinder(nuc_data, build_dir):
    """Controller function for adding cinder data."""
    with tb.openFile(nuc_data, 'a') as f:
        if hasattr(f.root, 'neutron') and hasattr(f.root.neutron, 'cinder_xs') and hasattr(f.root.neutron, 'cinder_fission_products'):
            return

    # First grab the atomic abundance data
    grab_cinder_dat(build_dir)

    # Add energy groups to file
    print "Adding cinder data..."
    print "  energy group boundaries."
    make_mg_group_structure(nuc_data, build_dir)

    # Add neutron absorption to file
    print "  neutron absorption cross sections."
    make_mg_absorption(nuc_data, build_dir)

    # Add fission to file
    print "  neutron fission cross sections."
    make_mg_fission(nuc_data, build_dir)

    # Add gamma decay spectrum to file
    print "  gamma decay spectra."
    make_mg_gamma_decay(nuc_data, build_dir)

    # Add neutron info table
    print "  neutron fission product info."
    make_neutron_fp_info(nuc_data, build_dir)

    # Add neutron yield table
    print "  neutron fission product yields."
    make_neutron_fp_yields(nuc_data, build_dir)

    # Add photon info table
    print "  photofission product info."
    make_photon_fp_info(nuc_data, build_dir)

    # Add neutron yield table
    print "  photofission product yields."
    make_photon_fp_yields(nuc_data, build_dir)

    print "...finished with cinder."
