#!/usr/bin/env python
import argparse
import yaml
import io

import numpy as np
from pyne.mesh import Mesh
from pyne.partisn import write_partisn_input, isotropic_vol_source
from pyne.dagmc import discretize_geom, load
from pyne import nucname
from pyne.bins import pointwise_collapse


config_filename = 'config.yml'

config = \
    """
# Optional step to assess all materials in geometry for compatibility with
# SNILB criteria
step0:

# Prepare PARTISN input file for adjoint photon transport
step1:
    # Path to hdf5 geometry file for photon transport
    geom_file:
    # Volume ID of adjoint photon source cell on
    # DAGMC input [Trelis/Cubit .sat file]
    src_cell:
    # Volume [cm^3] of source cell (detector)
    src_vol:
    # Define uniformly spaced, rectagular mesh that covers entire geometry:
    # Define origin of the mesh (three entries, one per spatial dimension)
    # Define locations of the coarse meshes in each direction
    # and the number of fine meshes within corresponding coarse meshes
    # Supported: Only one entry per _mesh and _ints for a uniformly 
    # spaced mesh
    # Separate values with blank space.
    origin:
    xmesh:
    xints:
    ymesh:
    yints:
    zmesh:
    zints:

# Calculate T matrix for each material
step2:

# Calculate adjoint neutron source
step3:

# Prepare PARTISN input for adjoint neutron transport
step4:

# Generate Monte Carlo variance reduction parameters
# (biased source and weight windows)
step5:


"""


def setup():
    """ This function generates a blank config.yml file for the user to 
    fill in with problem specific values. 
    """
    with open(config_filename, 'w') as f:
        f.write(config)
    print('File "{}" has been written'.format(config_filename))
    print('Fill out the fields in this file then run ">> gtcadis.py step1"')


def _names_dict():
    names = {'h1': 'h1', 'h2': 'd', 'h3': 'h3', 'he3': 'he3',
             'he4': 'he4', 'li6': 'li6', 'li7': 'li7',
             'be9': 'be9', 'b10': 'b10', 'b11': 'b11',
             'c12': 'c12', 'n14': 'n14', 'n15': 'n15',
             'o16': 'o16', 'f19': 'f19', 'na23': 'na23',
             'mg': 'mgnat', 'al27': 'al27', 'si28': 'si28',
             'si29': 'si29', 'si30': 'si30', 'p31': 'p31',
             's': 'snat', 'cl35': 'cl35', 'cl37': 'cl37',
             'k': 'knat', 'ca': 'canat', 'ti46': 'ti46',
             'ti47': 'ti47', 'ti48': 'ti48', 'ti49': 'ti49',
             'ti50': 'ti50', 'v': 'vnat', 'cr50': 'cr50',
             'cr52': 'cr52', 'cr53': 'cr53', 'cr54': 'cr54',
             'mn55': 'mn55', 'fe54': 'fe54', 'fe56': 'fe56',
             'fe57': 'fe57', 'fe58': 'fe58', 'co59': 'co59',
             'ni58': 'ni58', 'ni60': 'ni60', 'ni61': 'ni61',
             'ni62': 'ni62', 'ni64': 'ni64', 'cu63': 'cu63',
             'cu65': 'cu65', 'ga': 'ganat', 'zr': 'zrnat',
             'nb93': 'nb93', 'mo92': 'mo92', 'mo94': 'mo94',
             'mo95': 'mo95', 'mo96': 'mo96', 'mo97': 'mo97',
             'mo98': 'mo98', 'mo100': 'mo100', 'sn': 'snnat',
             'ta181': 'ta181', 'w182': 'w182', 'w183': 'w183',
             'w184': 'w184', 'w186': 'w186', 'au197': 'au197',
             'pb206': 'pb206', 'pb207': 'pb207', 'pb208': 'pb208',
             'bi209': 'bi209'}

    names_dict = {nucname.id(key): value for key, value in names.iteritems()}
    return names_dict


def _cards(source):
    cards = {"block1": {"isn": 16,
                        "maxscm": '3E8',
                        "maxlcm": '6E8',
                        },
             "block3": {"lib": "xsf21",  # name of cross section library
                        "lng": 175,
                        "maxord": 5,
                        "ihm": 227,
                        "iht": 10,
                        "ihs": 11,
                        "ifido": 1,
                        "ititl": 1,
                        "i2lp1": 0,
                        "savbxs": 0,
                        "kwikrd": 0
                        },
             "block5": {"source": source,
                        "ith": 1,
                        "isct": 5}
             }
    return cards


def step1(cfg):
    """ This function writes the PARTISN input file for the adjoint photon 
    transport   
    Parameters
    ----------
    cfg : dictionary
        User input for step 1 from the config.yml file
    """
    # Get user-input from config file
    geom = cfg['geom_file']
    cells = [cfg['src_cell']]
    src_vol = [float(cfg['src_vol'])]

    try:
        origin_x, origin_y, origin_z = cfg['origin'].split(' ')
    except:
        print("Too few entries in origin location")

    xmesh = cfg['xmesh']
    xints = cfg['xints']
    ymesh = cfg['ymesh']
    yints = cfg['yints']
    zmesh = cfg['zmesh']
    zints = cfg['zints']

    # Create structured mesh
    sc = [np.linspace(float(origin_x), float(xmesh), float(xints) + 1),
          np.linspace(float(origin_y), float(ymesh), float(yints) + 1),
          np.linspace(float(origin_z), float(zmesh), float(zints) + 1)]
    m = Mesh(structured=True, structured_coords=sc)
    m.write_hdf5("blank_mesh.h5m")

    # Generate 42 photon energy bins [eV]
    #  First bin has been replaced with 1 for log interpolation
    photon_bins = np.array([1e-6, 0.01, 0.02, 0.03, 0.045, 0.06, 0.07, 0.075, 0.1, 0.15,
                            0.2, 0.3, 0.4, 0.45, 0.51, 0.512, 0.6, 0.7, 0.8, 1, 1.33, 1.34,
                            1.5, 1.66, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 10,
                            12, 14, 20, 30, 50])
    # ICRP 74 flux-to-dose conversion factors in pico-Sv/s per photon flux
    de = np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.15, 0.2, 0.3,
                   0.4, 0.5, 0.6, 0.8, 1, 2, 4, 6, 8, 10])
    df = np.array([0.0485, 0.1254, 0.205, 0.2999, 0.3381, 0.3572, 0.378, 0.4066, 0.4399, 0.5172,
                   0.7523, 1.0041, 1.5083, 1.9958, 2.4657, 2.9082, 3.7269, 4.4834, 7.4896,
                   12.0153, 15.9873, 19.9191, 23.76])
    # Convert to Sv/s per photon FLUX
    pico = 1.0e-12
    df = df * pico
    # Convert pointwise data to group data for log interpolation
    photon_spectrum = pointwise_collapse(
        photon_bins, de, df, logx=True, logy=True)
    #  Anything below 0.01 MeV should be assigned the DF value of 0.01 MeV
    photon_spectrum[0] = df[0]
    # Total number of groups is 217 (42 photon + 175 neutron)
    spectra = [np.append(photon_spectrum, np.zeros(175))]
    # The spectrum is normalized by PyNE, so we need to mutliply by the sum of
    # intensities in the spectrum.
    # Additionally, we divide by the volume of the source cell in order to get
    # source density.
    intensities = [np.sum(spectra) / src_vol]

    # Load geometry into DAGMC
    load(geom)
    # Generate isotropic photon volume source
    source, dg = isotropic_vol_source(geom, m, cells, spectra, intensities)

    # PARTISN input
    ngroup = 217  # total number of energy groups
    cards = _cards(source)  # block 1, 3, 5 input values
    names_dict = _names_dict()  # dictionary of isotopes (PyNE nucids to bxslib names)

    write_partisn_input(
        m,
        geom,
        ngroup,
        cards=cards,
        dg=dg,
        names_dict=names_dict,
        data_hdf5path="/materials",
        nuc_hdf5path="/nucid",
        fine_per_coarse=1)


def main():
    """ This function manages the setup and steps 1-5 for the GT-CADIS workflow.
    """

    gtcadis_help = ('This script automates the GT-CADIS process of \n'
                    'producing variance reduction parameters to optimize the\n'
                    'neutron transport step of the Rigorous 2-Step (R2S) method.\n')
    setup_help = ('Prints the file "config.yml" to be\n'
                  'filled in by the user.\n')
    step1_help = 'Creates the PARTISN input file for adjoint photon transport.'
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help=gtcadis_help, dest='command')

    setup_parser = subparsers.add_parser('setup', help=setup_help)
    step1_parser = subparsers.add_parser('step1', help=step1_help)

    args, other = parser.parse_known_args()
    if args.command == 'setup':
        setup()
    else:
        with open(config_filename, 'r') as f:
            cfg = yaml.load(f)

    if args.command == 'step1':
        step1(cfg['step1'])


if __name__ == '__main__':
    main()
